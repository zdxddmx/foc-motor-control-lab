/* main.c — FOC 控制环仿真：启动 -> 加速 -> 加负载 -> 稳态
 *
 *   ./focsim pmsm     表贴式 PMSM，id=0 控制
 *   ./focsim bldc     凸极 BLDC 参数对比
 *   ./focsim ff       对比"有/无前馈解耦"的电流环跟随误差
 *   ./focsim fw       对比"有/无弱磁"的高速能力
 */
#include "foc.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static const char *scen = "pmsm";

static void header(const char *title)
{
    printf("\n===== %s =====\n", title);
}

/* 跑一段仿真，返回统计。speed_ref 为机械角速度 rad/s，负载在 t_load 时加入 */
static FocStats run(const MotorParams *p, double speed_ref, double t_load,
                    double tl, double total_t, int use_decouple, int use_fw,
                    double omega_base, double *speed_end, double *iq_end)
{
    FocConfig cfg = foc_default_config();
    FocCtrl c;
    FocStats st;
    MotorState m;
    SensorState s;
    InverterState inv;

    foc_init(&c, &cfg);
    motor_init(&m);
    sensor_init(&s, cfg.speed_lpf_alpha);
    inverter_init(&inv, 1.0);          /* 1 us 死区 */
    foc_stats_reset(&st);

    const double dt = cfg.ts;          /* 20 kHz */
    const int n = (int)(total_t / dt);
    int32_t enc = 0;
    c.speed_ref = speed_ref;
    m.tl = 0.0;

    for (int i = 0; i < n; i++) {
        const double t = i * dt;
        if (t >= t_load) m.tl = tl;

        /* 弱磁（可选） */
        c.i_ref.d = use_fw ? foc_field_weakening(m.omega_e, omega_base, 8.0) : 0.0;

        /* 编码器计数按机械角推进（4000 counts/rev） */
        enc = (int32_t)(m.theta_e / (FOC_TWO_PI) * 4000.0 * p->pole_pairs);
        /* 用真实机械角速度产生计数增长，模拟编码器 */
        static double theta_m_acc;
        if (i == 0) theta_m_acc = 0.0;
        theta_m_acc += m.omega_m * dt;
        enc = (int32_t)(theta_m_acc / FOC_TWO_PI * 4000.0);
        sensor_update_from_encoder(&s, enc, dt, p->pole_pairs);

        /* 不做解耦时：把前馈量清零（用于对比实验） */
        Dq save_u = c.u_ref;
        foc_step(&c, &cfg, p, &m, &s, &inv, dt);
        if (!use_decouple) {
            /* 撤掉解耦项：重新用纯 PI 输出驱动 */
            c.u_ref.d = c.u_ref.d - (-m.omega_e * p->lq * c.i_fb.q);
            c.u_ref.q = c.u_ref.q - ( m.omega_e * (p->ld * c.i_fb.d + p->flux));
            svpwm_calc(c.u_ref, m.theta_e, cfg.udc, &c.pwm);
            inverter_apply(&inv, &c.pwm, cfg.udc);
            const Abc uabc = { inv.ua, inv.ub, inv.uc };
            const Dq udq = park(clarke(uabc), m.theta_e);
            m.ud = udq.d;
            m.uq = udq.q;
        }
        (void)save_u;

        motor_step(&m, p, dt);
        foc_stats_update(&st, &m, p, &inv, cfg.udc);
    }

    if (speed_end) *speed_end = m.omega_m;
    if (iq_end) *iq_end = m.iq;
    return st;
}

static void report(const char *tag, const FocStats *st, double speed_ref,
                   double speed_end, double iq_end)
{
    printf("  %-14s 转速 %8.2f / 目标 %8.2f rad/s（稳态误差 %+6.2f%%）\n",
           tag, speed_end, speed_ref,
           speed_ref != 0.0 ? (speed_end - speed_ref) * 100.0 / speed_ref : 0.0);
    printf("  %-14s iq均值=%7.3f A  id均值=%7.3f A  转矩均值=%7.4f N·m（脉动 %.4f）\n",
           "", st->iq_mean, st->id_mean, st->torque_mean, st->torque_ripple);
    printf("  %-14s 相电流 RMS A/B/C = %.3f / %.3f / %.3f A\n",
           "", sqrt(st->ia_rms), sqrt(st->ib_rms), sqrt(st->ic_rms));
    printf("  %-14s 效率≈%.1f%%（dq 轴口径）  桥臂直通=%u  仿真步数=%llu\n",
           "", st->efficiency * 100.0, st->shoot_through,
           (unsigned long long)st->samples);
}

int main(int argc, char **argv)
{
    if (argc > 1) scen = argv[1];
    const double T = 0.30;                 /* 每段仿真 0.3 s */

    if (strcmp(scen, "pmsm") == 0) {
        header("表贴式 PMSM · id=0 控制（启动 → 加速 → 加负载）");
        const MotorParams p = motor_default_pmsm();
        double se = 0, iq = 0;
        FocStats st = run(&p, 300.0, 0.12, 0.05, T, 1, 0, 1e9, &se, &iq);
        report("PMSM", &st, 300.0, se, iq);
        printf("  参数: Rs=%.2fΩ Ld=Lq=%.4fH ψf=%.4fWb P=%.0f  Udc=48V\n",
               p.rs, p.ld, p.flux, p.pole_pairs);
        printf("  电流环 20 kHz，死区 1 us，编码器 1000 线四倍频\n");
        return 0;
    }

    if (strcmp(scen, "bldc") == 0) {
        header("凸极 BLDC · 参数对比（Ld≠Lq）");
        const MotorParams p = motor_default_bldc();
        double se = 0, iq = 0;
        FocStats st = run(&p, 300.0, 0.12, 0.05, T, 1, 0, 1e9, &se, &iq);
        report("BLDC", &st, 300.0, se, iq);
        printf("  参数: Rs=%.2fΩ Ld=%.4f Lq=%.4f ψf=%.4f P=%.0f\n",
               p.rs, p.ld, p.lq, p.flux, p.pole_pairs);
        /* 六步换相表验证 */
        SixStepOut ss;
        printf("\n  六步换相表（霍尔 1~6）：\n");
        for (uint8_t h = 1; h <= 6; h++) {
            bldc_six_step(h, 0.5, &ss);
            printf("    hall=%u -> 上桥 %c / 下桥 %c / 悬空 %c  下一步=%u\n",
                   h, 'A' + ss.phase_high, 'A' + ss.phase_low,
                   'A' + ss.floating, bldc_next_step(h));
        }
        return 0;
    }

    if (strcmp(scen, "ff") == 0) {
        header("前馈解耦对比（高速下电流环跟随误差）");
        const MotorParams p = motor_default_pmsm();
        const double wref = 800.0;                 /* 高速 */
        double se1 = 0, iq1 = 0, se2 = 0, iq2 = 0;
        FocStats a = run(&p, wref, 1e9, 0.0, T, 1, 0, 1e9, &se1, &iq1);
        FocStats b = run(&p, wref, 1e9, 0.0, T, 0, 0, 1e9, &se2, &iq2);
        report("有解耦", &a, wref, se1, iq1);
        report("无解耦", &b, wref, se2, iq2);
        printf("\n  结论：高转速下不做前馈解耦时，q 轴电流被反电动势与交叉耦合项拖住，\n");
        printf("        转速稳态误差与电流跟随误差都会明显变大。\n");
        return 0;
    }

    if (strcmp(scen, "fw") == 0) {
        header("弱磁扩速对比");
        const MotorParams p = motor_default_pmsm();
        const double wbase = 320.0;
        const double wref = 460.0;                 /* 超过基速 */
        double se1 = 0, iq1 = 0, se2 = 0, iq2 = 0;
        FocStats a = run(&p, wref, 1e9, 0.0, T, 1, 1, wbase, &se1, &iq1);
        FocStats b = run(&p, wref, 1e9, 0.0, T, 1, 0, wbase, &se2, &iq2);
        report("有弱磁", &a, wref, se1, iq1);
        report("无弱磁", &b, wref, se2, iq2);
        printf("\n  结论：超过基速后母线电压不足以维持电流环，弱磁给负 id 削弱磁链，\n");
        printf("        才能继续升速（代价是同样的转矩需要更大的电流）。\n");
        return 0;
    }

    fprintf(stderr, "未知场景: %s（可选 pmsm / bldc / ff / fw）\n", scen);
    return 2;
}
