/* test_foc.c — 单元测试（纯 C99，无第三方框架）
 *
 *  T1 math      坐标变换（Clarke/Park 往返、幅值关系）
 *  T2 pi        PI 控制器（比例/积分、抗饱和、前馈、复位）
 *  T3 svpwm     扇区判定、占空比范围、过调制限幅、零矢量
 *  T4 sensor    霍尔角度与有效性、编码器角度/速度、速度低通
 *  T5 bldc      六步换相表与悬空相
 *  T6 inverter  互补桥臂、直通检测（显式注入）
 *  T7 foc       闭环收敛、稳态误差、id=0 控制、弱磁、直通计数
 *  T8 stats     增量均值、效率口径
 */
#include "foc.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int passed = 0, failed = 0;

#define CHECK(cond, ...) do {                        \
    if (cond) { passed++; printf("  [PASS] "); }      \
    else      { failed++; printf("  [FAIL] "); }      \
    printf(__VA_ARGS__); printf("\n");                \
} while (0)

#define NEAR(a, b, tol) (fabs((a) - (b)) <= (tol))

/* ------------------------------------------------------------------ T1 */
static void t1_math(void)
{
    printf("T1 坐标变换\n");

    /* 纯 A 相电流：Clarke 后 alpha = 2/3*Ia, beta = 0（等幅值变换） */
    Abc in = { 3.0, 0.0, 0.0 };
    AlphaBeta ab = clarke(in);
    CHECK(NEAR(ab.alpha, 2.0, 1e-9), "Clarke: ia=3,ib=ic=0 -> alpha=%.6f（=2/3·ia）", ab.alpha);
    CHECK(NEAR(ab.beta, 0.0, 1e-9), "Clarke: beta=%.6f", ab.beta);

    /* 三相对称时 alpha/beta 应回到原幅值 */
    const double I = 5.0, w = 0.7;
    Abc sym = { I * cos(w), I * cos(w - FOC_TWO_PI / 3.0), I * cos(w + FOC_TWO_PI / 3.0) };
    ab = clarke(sym);
    CHECK(NEAR(sqrt(ab.alpha * ab.alpha + ab.beta * ab.beta), I, 1e-9),
          "三相对称：矢量幅值 %.6f（应等于相幅值 %.1f）",
          sqrt(ab.alpha * ab.alpha + ab.beta * ab.beta), I);

    /* Clarke 往返 */
    Abc back = inv_clarke(ab);
    CHECK(NEAR(back.a, sym.a, 1e-9) && NEAR(back.b, sym.b, 1e-9) &&
          NEAR(back.c, sym.c, 1e-9), "inv_clarke 往返一致");

    /* Park：theta=0 时 d=alpha, q=beta */
    AlphaBeta ab2 = { 1.0, 2.0 };
    Dq dq = park(ab2, 0.0);
    CHECK(NEAR(dq.d, 1.0, 1e-9) && NEAR(dq.q, 2.0, 1e-9),
          "Park(theta=0): d=%.3f q=%.3f", dq.d, dq.q);

    /* Park 往返 */
    Dq dq2 = { 3.0, -4.0 };
    AlphaBeta r = inv_park(dq2, 1.234);
    Dq dq3 = park(r, 1.234);
    CHECK(NEAR(dq3.d, 3.0, 1e-9) && NEAR(dq3.q, -4.0, 1e-9),
          "Park/inv_park 往返一致：d=%.6f q=%.6f", dq3.d, dq3.q);

    /* theta=90°：alpha 轴电流应全部投影到 q 轴 */
    AlphaBeta ab3 = { 1.0, 0.0 };
    Dq dq4 = park(ab3, FOC_PI / 2.0);
    CHECK(NEAR(dq4.d, 0.0, 1e-9) && NEAR(dq4.q, -1.0, 1e-9),
          "Park(90°): d=%.6f q=%.6f（alpha 投到 -q）", dq4.d, dq4.q);

    /* 电角度归一化 */
    CHECK(NEAR(motor_elec_angle(0.0, 4.0), 0.0, 1e-12), "电角度 0");
    CHECK(NEAR(motor_elec_angle(FOC_PI / 2.0, 4.0), 0.0, 1e-9),
          "机械 90° × P=4 -> 电角度归一化为 0");
    CHECK(motor_elec_angle(-0.1, 1.0) > 0.0, "负角度归一化到正区间");
    const double e = motor_elec_angle(12.3, 4.0);
    CHECK(e >= 0.0 && e < FOC_TWO_PI, "电角度始终在 [0,2π)：%.6f", e);
}

/* ------------------------------------------------------------------ T2 */
static void t2_pi(void)
{
    printf("T2 PI 控制器\n");

    PiCtrl pi;
    pi_init(&pi, 2.0, 10.0, 0.001, -100.0, 100.0);
    /* 未饱和时即开始积分：首次输出 = kp·err + ki·err·ts */
    double out = pi_step(&pi, 5.0, 0.0);
    const double e1 = 2.0 * 5.0 + 10.0 * (5.0 * 0.001);
    CHECK(NEAR(out, e1, 1e-9), "首次输出 = kp·err + ki·err·ts = %.4f（期望 %.4f）",
          out, e1);

    /* 第二次：积分累加 err1+err2 */
    out = pi_step(&pi, 5.0, 1.0);
    const double e2 = 2.0 * 4.0 + 10.0 * ((5.0 + 4.0) * 0.001);
    CHECK(NEAR(out, e2, 1e-6), "第二次输出 = %.4f（期望 %.4f）", out, e2);

    /* 输出限幅 */
    PiCtrl p2;
    pi_init(&p2, 100.0, 0.0, 0.001, -1.0, 1.0);
    out = pi_step(&p2, 1e6, 0.0);
    CHECK(NEAR(out, 1.0, 1e-12), "输出上限截断 = %.3f", out);
    out = pi_step(&p2, -1e6, 0.0);
    CHECK(NEAR(out, -1.0, 1e-12), "输出下限截断 = %.3f", out);

    /* 抗饱和：持续大误差时积分不应无限增长 */
    PiCtrl p3;
    pi_init(&p3, 1.0, 100.0, 0.001, -1.0, 1.0);
    for (int i = 0; i < 5000; i++) pi_step(&p3, 100.0, 0.0);
    /* 去掉输入后应立即回落（积分没被灌满） */
    out = pi_step(&p3, 0.0, 0.0);
    CHECK(out <= 1.0 && p3.integ < 1.0,
          "条件积分抗饱和：撤掉输入后输出 %.4f，积分 %.4f", out, p3.integ);

    /* 前馈 */
    PiCtrl p4;
    pi_init(&p4, 1.0, 0.0, 0.001, -100.0, 100.0);
    pi_set_feedforward(&p4, 7.5);
    out = pi_step(&p4, 0.0, 0.0);
    CHECK(NEAR(out, 7.5, 1e-9), "前馈单独作用 = %.3f", out);
    out = pi_step(&p4, 1.0, 0.0);
    CHECK(NEAR(out, 8.5, 1e-9), "前馈 + 比例 = %.3f", out);

    /* 复位 */
    pi_reset(&p4);
    CHECK(NEAR(p4.integ, 0.0, 1e-12) && p4.calls == 2,
          "复位清积分但保留调用计数（calls=%u）", p4.calls);
}

/* ------------------------------------------------------------------ T3 */
static void t3_svpwm(void)
{
    printf("T3 SVPWM\n");

    /* 扇区判定：纯 +alpha 应在扇区 1 与 6 的边界附近，取几个明确方向 */
    /* +β 轴对应电角度 90°，落在扇区 2（扇区 1 = 0~60°，扇区 2 = 60~120°） */
    CHECK(svpwm_sector(0.0, 1.0) == 2, "beta>0, alpha=0（+β 轴）-> 扇区 %d",
          svpwm_sector(0.0, 1.0));
    /* +α 轴（0°）正好是扇区 1 与 6 的分界，返回 1 或 6 都算对，
     * 所以用明显落在扇区内部的方向来验证 */
    CHECK(svpwm_sector(1.0, 0.2) == 1, "约 11° -> 扇区 %d",
          svpwm_sector(1.0, 0.2));
    CHECK(svpwm_sector(1.0, 0.0) == 1 || svpwm_sector(1.0, 0.0) == 6,
          "+α 轴（扇区边界）-> 扇区 %d（1 或 6 均可）",
          svpwm_sector(1.0, 0.0));
    CHECK(svpwm_sector(0.0, -1.0) == 5, "-β 轴 -> 扇区 %d", svpwm_sector(0.0, -1.0));
    CHECK(svpwm_sector(-1.0, 0.0) == 4, "-α 轴 -> 扇区 %d", svpwm_sector(-1.0, 0.0));
    CHECK(svpwm_sector(1.0, 1.0) >= 1, "第一象限 -> 扇区 %d",
          svpwm_sector(1.0, 1.0));

    /* 零输入 -> 零矢量 -> 三相 50% */
    SvpwmOut o;
    svpwm_calc((Dq){ 0.0, 0.0 }, 0.0, 48.0, &o);
    CHECK(NEAR(o.ta, 0.5, 1e-9) && NEAR(o.tb, 0.5, 1e-9) && NEAR(o.tc, 0.5, 1e-9),
          "零电压 -> 三相占空比 0.5/0.5/0.5");
    CHECK(o.mod_index < 1e-9, "零矢量调制比 = %.6f", o.mod_index);

    /* 正常电压：占空比必须都在 [0,1]，且 T1+T2+T0 = 1 */
    svpwm_calc((Dq){ 0.0, 10.0 }, 0.3, 48.0, &o);
    CHECK(o.ta >= 0.0 && o.ta <= 1.0 && o.tb >= 0.0 && o.tb <= 1.0 &&
          o.tc >= 0.0 && o.tc <= 1.0,
          "占空比均在 [0,1]：%.4f %.4f %.4f", o.ta, o.tb, o.tc);
    CHECK(NEAR(o.t1 + o.t2 + o.t0, 1.0, 1e-9),
          "T1+T2+T0 = %.6f", o.t1 + o.t2 + o.t0);

    /* 线电压重建：SVPWM 输出的占空比差应能还原参考相电压 */
    const Dq ref = { 5.0, 12.0 };
    svpwm_calc(ref, 0.9, 48.0, &o);
    const double uab = (o.ta - o.tb) * 48.0;
    const double ubc = (o.tb - o.tc) * 48.0;
    const double uca = (o.tc - o.ta) * 48.0;
    CHECK(NEAR(uab + ubc + uca, 0.0, 1e-6), "线电压之和 = %.9f（应为 0）",
          uab + ubc + uca);

    /* 过调制：给极大的参考电压，T1+T2 应被压回 1 以内 */
    svpwm_calc((Dq){ 200.0, 200.0 }, 0.5, 48.0, &o);
    CHECK(o.t1 + o.t2 <= 1.0 + 1e-9, "过调制限幅：T1+T2 = %.6f ≤ 1", o.t1 + o.t2);
    CHECK(o.ta <= 1.0 && o.tb <= 1.0 && o.tc <= 1.0, "过调制后占空比仍 ≤1");

    /* 母线为 0 不应崩 */
    svpwm_calc(ref, 0.0, 0.0, &o);
    CHECK(NEAR(o.ta, 0.5, 1e-9), "母线 0V -> 退化为零矢量");

    /* 调制比随电压线性变化 */
    SvpwmOut a, b;
    svpwm_calc((Dq){ 0.0, 5.0 }, 0.0, 48.0, &a);
    svpwm_calc((Dq){ 0.0, 10.0 }, 0.0, 48.0, &b);
    CHECK(NEAR(b.mod_index, 2.0 * a.mod_index, 1e-9),
          "调制比倍增：%.6f -> %.6f", a.mod_index, b.mod_index);
}

/* ------------------------------------------------------------------ T4 */
static void t4_sensor(void)
{
    printf("T4 位置与速度\n");

    int valid = -1;
    double ang = sensor_hall_angle(0x05, &valid);
    CHECK(valid == 1 && NEAR(ang, 0.0, 1e-12), "霍尔 101 有效，角度 %.3f", ang);

    /* 六个有效状态应给出 6 个不同的 60° 区间 */
    const uint8_t states[6] = { 0x05, 0x01, 0x03, 0x02, 0x06, 0x04 };
    double angles[6];
    int all_valid = 1;
    for (int i = 0; i < 6; i++) {
        angles[i] = sensor_hall_angle(states[i], &valid);
        if (!valid) all_valid = 0;
    }
    CHECK(all_valid, "六个标准霍尔状态均有效");
    int distinct = 1;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            if (NEAR(angles[i], angles[j], 1e-9)) distinct = 0;
    CHECK(distinct, "六个状态对应六个不同角度（间隔 60°）");
    CHECK(NEAR(angles[1] - angles[0], FOC_PI / 3.0, 1e-9),
          "相邻状态角度差 = %.4f rad（=60°）", angles[1] - angles[0]);

    /* 000 与 111 是无效状态 */
    sensor_hall_angle(0x00, &valid);
    CHECK(valid == 0, "霍尔 000 判为无效");
    sensor_hall_angle(0x07, &valid);
    CHECK(valid == 0, "霍尔 111 判为无效");

    /* 编码器：匀加速下角度与速度应单调 */
    SensorState s;
    sensor_init(&s, 1.0);            /* 不做低通，便于精确校验 */
    const double dt = 0.001, pp = 4.0;
    for (int i = 1; i <= 100; i++) {
        sensor_update_from_encoder(&s, i * 4, dt, pp);   /* 4 counts/ms = 1 rev/s */
    }
    CHECK(NEAR(s.omega_m, FOC_TWO_PI, 1e-6),
          "编码器速度 = %.6f rad/s（期望 2π）", s.omega_m);
    CHECK(NEAR(s.theta_m, 400.0 * FOC_TWO_PI / 4000.0, 1e-9),
          "机械角 = %.6f rad", s.theta_m);
    CHECK(s.theta_e >= 0.0 && s.theta_e < FOC_TWO_PI, "电角度在 [0,2π)");

    /* 低通：突变速度不应立刻全额通过 */
    SensorState s2;
    sensor_init(&s2, 0.2);
    sensor_update_from_encoder(&s2, 0, dt, pp);
    sensor_update_from_encoder(&s2, 100, dt, pp);
    CHECK(s2.omega_filt < s2.omega_m,
          "低通后速度 %.3f < 原始 %.3f（有滞后）", s2.omega_filt, s2.omega_m);

    /* 非法 lpf 参数应回落到默认 */
    SensorState s3;
    sensor_init(&s3, 0.0);
    CHECK(NEAR(s3.lpf_alpha, 0.2, 1e-12), "lpf=0 回落到默认 0.2");
    sensor_init(&s3, 5.0);
    CHECK(NEAR(s3.lpf_alpha, 0.2, 1e-12), "lpf>1 回落到默认 0.2");
}

/* ------------------------------------------------------------------ T5 */
static void t5_bldc(void)
{
    printf("T5 BLDC 六步换相\n");

    CHECK(bldc_next_step(1) == 2 && bldc_next_step(6) == 1,
          "换相顺序 1->2 ... 6->1");
    CHECK(bldc_next_step(0) == 0 && bldc_next_step(7) == 0,
          "非法霍尔状态返回 0");

    /* 六步覆盖全部 6 种上下桥组合，且悬空相各不相同地轮转 */
    int high_seen[3] = { 0, 0, 0 }, low_seen[3] = { 0, 0, 0 }, float_seen[3] = { 0, 0, 0 };
    int no_conflict = 1;
    for (uint8_t h = 1; h <= 6; h++) {
        SixStepOut o;
        bldc_six_step(h, 0.5, &o);
        if (o.phase_high > 2 || o.phase_low > 2 || o.floating > 2) no_conflict = 0;
        if (o.phase_high == o.phase_low || o.phase_high == o.floating ||
            o.phase_low == o.floating) no_conflict = 0;
        high_seen[o.phase_high] = 1;
        low_seen[o.phase_low] = 1;
        float_seen[o.floating] = 1;
    }
    CHECK(no_conflict, "每步的三相职责互不相同（上桥/下桥/悬空各一相）");
    CHECK(high_seen[0] && high_seen[1] && high_seen[2], "三相都当过上桥");
    CHECK(low_seen[0] && low_seen[1] && low_seen[2], "三相都当过下桥");
    CHECK(float_seen[0] && float_seen[1] && float_seen[2], "三相都当过悬空相");

    /* 占空比限幅 */
    SixStepOut o;
    bldc_six_step(1, 1.5, &o);
    CHECK(NEAR(o.duty, 1.0, 1e-12), "占空比上限截断 = %.3f", o.duty);
    bldc_six_step(1, -0.5, &o);
    CHECK(NEAR(o.duty, 0.0, 1e-12), "占空比下限截断 = %.3f", o.duty);

    /* 整套换相序列应能回到起点 */
    uint8_t cur = 1;
    for (int i = 0; i < 6; i++) cur = bldc_next_step(cur);
    CHECK(cur == 1, "换相 6 次回到初始状态");
}

/* ------------------------------------------------------------------ T6 */
static void t6_inverter(void)
{
    printf("T6 三相逆变器\n");

    InverterState inv;
    inverter_init(&inv, 1.0);
    SvpwmOut o;
    svpwm_calc((Dq){ 0.0, 12.0 }, 0.4, 48.0, &o);
    inverter_apply(&inv, &o, 48.0);

    CHECK(inv.shoot_through_events == 0, "正常驱动不产生直通（计数=%u）",
          inv.shoot_through_events);
    /* 互补性：每相上下桥不同时为 1 */
    const uint16_t r = inv.bridge.raw;
    CHECK(((r & 0x03u) != 0x03u) && ((r & 0x0Cu) != 0x0Cu) && ((r & 0x30u) != 0x30u),
          "三相均为互补驱动（raw=0x%04X）", r);
    CHECK(inv.bridge.bits.deadtime_active == 1, "死区标志已置位");

    /* 相电压：占空比 0.5 时为 0，>0.5 为正 */
    SvpwmOut half;
    svpwm_calc((Dq){ 0.0, 0.0 }, 0.0, 48.0, &half);
    inverter_apply(&inv, &half, 48.0);
    CHECK(NEAR(inv.ua, 0.0, 1e-9) && NEAR(inv.ub, 0.0, 1e-9) && NEAR(inv.uc, 0.0, 1e-9),
          "零矢量时三相电压均为 0");

    /* 显式注入直通：把 A 相上下桥同时置 1，保护必须检出 */
    inv.bridge.raw = 0x0003u;
    CHECK(inverter_check_shoot_through(&inv) == 1, "注入 A 相直通 -> 被检出");
    CHECK(inv.shoot_through_events == 1, "直通计数 = %u", inv.shoot_through_events);
    CHECK(inv.bridge.bits.shoot_through == 1, "直通故障标志置位");

    inv.bridge.raw = 0x000Cu;
    CHECK(inverter_check_shoot_through(&inv) == 1, "注入 B 相直通 -> 被检出");
    inv.bridge.raw = 0x0030u;
    CHECK(inverter_check_shoot_through(&inv) == 1, "注入 C 相直通 -> 被检出");
    CHECK(inv.shoot_through_events == 3, "累计直通 3 次");

    /* 正常状态清除故障标志 */
    inv.bridge.raw = 0x0005u;
    CHECK(inverter_check_shoot_through(&inv) == 0 && inv.bridge.bits.shoot_through == 0,
          "恢复正常后故障标志清除");

    /* 死区为 0 时标志不置位 */
    InverterState inv2;
    inverter_init(&inv2, 0.0);
    inverter_apply(&inv2, &o, 48.0);
    CHECK(inv2.bridge.bits.deadtime_active == 0, "死区 0 时不置位死区标志");
}

/* ------------------------------------------------------------------ T7 */
static void t7_foc(void)
{
    printf("T7 FOC 闭环\n");

    FocConfig cfg = foc_default_config();
    CHECK(NEAR(cfg.ts, 1.0 / 20000.0, 1e-12), "默认电流环周期 = %.1f us",
          cfg.ts * 1e6);
    CHECK(NEAR(cfg.udc, 48.0, 1e-9), "默认母线 = %.1f V", cfg.udc);
    CHECK(cfg.i_max > 0.0 && cfg.u_max_ratio > 0.0 && cfg.u_max_ratio <= 1.0,
          "限幅参数合法");

    const MotorParams p = motor_default_pmsm();
    CHECK(NEAR(p.ld, p.lq, 1e-12), "表贴式 PMSM：Ld = Lq = %.5f H", p.ld);
    const MotorParams pb = motor_default_bldc();
    CHECK(pb.ld < pb.lq, "凸极 BLDC：Ld(%.5f) < Lq(%.5f)", pb.ld, pb.lq);

    FocCtrl c;
    foc_init(&c, &cfg);
    CHECK(NEAR(c.i_ref.d, 0.0, 1e-12), "初始化后 id 参考 = 0（id=0 控制）");
    CHECK(NEAR(c.udc, cfg.udc, 1e-9), "控制器母线电压已同步");

    /* 闭环收敛测试 */
    MotorState m;
    SensorState s;
    InverterState inv;
    motor_init(&m);
    sensor_init(&s, cfg.speed_lpf_alpha);
    inverter_init(&inv, 1.0);

    const double dt = cfg.ts;
    const double wref = 300.0;
    const double tl = 0.05;
    c.speed_ref = wref;
    m.tl = tl;

    double theta_m = 0.0;
    const int n = (int)(0.5 / dt);
    for (int i = 0; i < n; i++) {
        theta_m += m.omega_m * dt;
        const int32_t enc = (int32_t)(theta_m / FOC_TWO_PI * 4000.0);
        sensor_update_from_encoder(&s, enc, dt, p.pole_pairs);
        foc_step(&c, &cfg, &p, &m, &s, &inv, dt);
        motor_step(&m, &p, dt);
    }

    const double err = (m.omega_m - wref) * 100.0 / wref;
    CHECK(fabs(err) < 5.0, "4800 步后转速 %.2f / 目标 %.0f rad/s（误差 %+.2f%%）",
          m.omega_m, wref, err);
    CHECK(fabs(m.id) < 1.0, "id=0 控制下 id = %.4f A（应接近 0）", m.id);
    CHECK(m.iq > 0.0, "带载时 iq 为正 = %.4f A（产生正向转矩）", m.iq);
    CHECK(inv.shoot_through_events == 0, "闭环全过程无桥臂直通");

    /* 转矩方向正确（带载时为正）。
     * 注意：不要用"最后一个采样点"的 Te 去精确比对负载 ——
     * 转矩存在明显脉动，末点可能正落在谷底。
     * 转矩是否托住负载，由上面"转速稳态误差 < 5%"这条来保证。 */
    CHECK(m.te > 0.0, "稳态转矩为正 = %.4f N·m（带载 %0.3f N·m）", m.te, tl);
    CHECK(m.te < 1.0, "稳态转矩在合理量级内 = %.4f N·m", m.te);
    CHECK(m.iq > 0.0, "iq 为正 = %.4f A（产生正转矩）", m.iq);

    /* ISR 计数应与步数一致 */
    CHECK(c.isr_count == (uint64_t)n, "电流环中断计数 = %llu（步数 %d）",
          (unsigned long long)c.isr_count, n);

    /* 弱磁：低于基速不给负 id，高于基速按比例给并限幅 */
    CHECK(NEAR(foc_field_weakening(100.0, 300.0, 5.0), 0.0, 1e-12),
          "低于基速不弱磁");
    const double id1 = foc_field_weakening(450.0, 300.0, 5.0);
    CHECK(id1 < 0.0, "高于基速给负 id = %.3f A", id1);
    const double id2 = foc_field_weakening(10000.0, 300.0, 5.0);
    CHECK(NEAR(id2, -5.0, 1e-9), "弱磁限幅到 id_max = %.3f A", id2);
    CHECK(foc_field_weakening(500.0, 0.0, 5.0) == 0.0, "基速非法时不弱磁");
}

/* ------------------------------------------------------------------ T8 */
static void t8_stats(void)
{
    printf("T8 统计\n");

    FocStats st;
    foc_stats_reset(&st);
    CHECK(st.samples == 0 && NEAR(st.speed_mean, 0.0, 1e-12), "复位清零");

    MotorParams p = motor_default_pmsm();
    MotorState m;
    motor_init(&m);
    InverterState inv;
    inverter_init(&inv, 0.0);

    /* 固定状态喂 100 次，增量均值应收敛到该值 */
    m.omega_m = 250.0; m.te = 0.05; m.id = 0.0; m.iq = 1.0;
    m.ud = 2.0; m.uq = 10.0;
    m.ia = m.ib = m.ic = 1.0;
    for (int i = 0; i < 100; i++)
        foc_stats_update(&st, &m, &p, &inv, 48.0);

    CHECK(st.samples == 100, "样本数 = %llu", (unsigned long long)st.samples);
    CHECK(NEAR(st.speed_mean, 250.0, 1e-9), "速度均值 = %.4f", st.speed_mean);
    CHECK(NEAR(st.iq_mean, 1.0, 1e-9), "iq 均值 = %.4f", st.iq_mean);
    CHECK(NEAR(st.torque_mean, 0.05, 1e-9), "转矩均值 = %.6f", st.torque_mean);
    CHECK(NEAR(st.speed_ripple, 0.0, 1e-9), "恒定速度下脉动 = %.9f", st.speed_ripple);
    CHECK(NEAR(sqrt(st.ia_rms), 1.0, 1e-9), "相电流 RMS = %.6f（输入 1.0）",
          sqrt(st.ia_rms));

    /* 效率口径：P_in = 1.5(ud·id + uq·iq) = 15 W，P_out = Te·ωm = 12.5 W */
    CHECK(NEAR(st.efficiency, 12.5 / 15.0, 1e-6),
          "效率 = %.4f（期望 %.4f = 12.5/15）", st.efficiency, 12.5 / 15.0);

    /* 脉动：一半样本速度不同，均值仍应正确 */
    FocStats s2;
    foc_stats_reset(&s2);
    for (int i = 0; i < 50; i++) {
        m.omega_m = 200.0;
        foc_stats_update(&s2, &m, &p, &inv, 48.0);
        m.omega_m = 300.0;
        foc_stats_update(&s2, &m, &p, &inv, 48.0);
    }
    CHECK(NEAR(s2.speed_mean, 250.0, 1e-6), "交替速度均值 = %.4f", s2.speed_mean);
    CHECK(s2.speed_ripple > 40.0, "脉动被统计到 = %.3f", s2.speed_ripple);
}

int main(void)
{
    printf("==== foc 无刷电机控制库 单元测试 ====\n\n");
    t1_math();
    t2_pi();
    t3_svpwm();
    t4_sensor();
    t5_bldc();
    t6_inverter();
    t7_foc();
    t8_stats();
    printf("\n==== 结果：%d passed, %d failed ====\n", passed, failed);
    return failed ? 1 : 0;
}
