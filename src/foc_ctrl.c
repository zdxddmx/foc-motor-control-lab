/* foc_ctrl.c — SVPWM、位置/速度估算、BLDC 六步换相、三相逆变器、FOC 主循环 */
#include "foc.h"
#include <math.h>
#include <string.h>

/* ================================================================== SVPWM
 * 空间矢量脉宽调制：把 dq 电压换成三相占空比。
 *   1) 逆 Park 得 (uα, uβ)
 *   2) 判断扇区（1~6）
 *   3) 算相邻两矢量作用时间 T1/T2 与零矢量 T0
 *   4) 七段式/五段式合成三相占空比
 * 相对 SPWM 的优点是母线电压利用率高约 15%。
 */
uint8_t svpwm_sector(double alpha, double beta)
{
    /* 用三个判断量定扇区，避免 atan2（车用 MCU 上更省时间） */
    const double a = beta;
    const double b = -0.5 * beta + FOC_SQRT3_2 * alpha;
    const double c = -0.5 * beta - FOC_SQRT3_2 * alpha;
    uint8_t s = 0;
    if (a > 0.0) s |= 0x01;
    if (b > 0.0) s |= 0x02;
    if (c > 0.0) s |= 0x04;
    switch (s) {
    case 0x03: return 1;      /* a>0, b>0 */
    case 0x01: return 2;      /* a>0, c<0 */
    case 0x05: return 3;      /* c<0, a>0 -> 由 b<0 区分 */
    case 0x04: return 4;
    case 0x06: return 5;
    case 0x02: return 6;
    default:   return 0;      /* 零矢量附近 */
    }
}

static double clamp01(double v)
{
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}

void svpwm_calc(Dq udq, double theta_e, double udc, SvpwmOut *out)
{
    memset(out, 0, sizeof(*out));
    const AlphaBeta u = inv_park(udq, theta_e);
    const double umag = sqrt(u.alpha * u.alpha + u.beta * u.beta);
    out->mod_index = (udc > 0.0) ? (umag / (udc * FOC_1_SQRT3)) : 0.0;

    const uint8_t sec = svpwm_sector(u.alpha, u.beta);
    out->sector = sec;
    if (sec == 0 || udc <= 0.0) {
        out->ta = out->tb = out->tc = 0.5;   /* 零矢量：三相 50% */
        return;
    }

    /* 基本矢量幅值 2/3·Udc 归一化后的作用时间 */
    const double k = FOC_SQRT3 / udc;        /* 归一化系数 */
    double t1 = 0.0, t2 = 0.0;
    switch (sec) {
    case 1: t1 = k * ( FOC_SQRT3_2 * u.alpha - 0.5 * u.beta);
            t2 = k * ( FOC_SQRT3_2 * u.alpha + 0.5 * u.beta); break;
    case 2: t1 = k * ( FOC_SQRT3_2 * u.alpha + 0.5 * u.beta);
            t2 = k * (-FOC_SQRT3_2 * u.alpha + 0.5 * u.beta); break;
    case 3: t1 = k * ( u.beta);
            t2 = k * (-u.alpha); break;
    case 4: t1 = k * (-FOC_SQRT3_2 * u.alpha + 0.5 * u.beta);
            t2 = k * (-FOC_SQRT3_2 * u.alpha - 0.5 * u.beta); break;
    case 5: t1 = k * (-FOC_SQRT3_2 * u.alpha - 0.5 * u.beta);
            t2 = k * ( FOC_SQRT3_2 * u.alpha - 0.5 * u.beta); break;
    case 6: t1 = k * (-u.beta);
            t2 = k * ( u.alpha); break;
    default: break;
    }
    if (t1 < 0.0) t1 = 0.0;
    if (t2 < 0.0) t2 = 0.0;
    /* 过调制保护：T1+T2 不得超过一个周期 */
    if (t1 + t2 > 1.0) {
        const double s = 1.0 / (t1 + t2);
        t1 *= s;
        t2 *= s;
    }
    const double t0 = 1.0 - t1 - t2;
    out->t1 = t1;
    out->t2 = t2;
    out->t0 = t0;

    /* 五段式：以零矢量为起点，按扇区把 T1/T2 分配到三相 */
    const double half0 = t0 * 0.5;
    double da = 0.0, db = 0.0, dc = 0.0;
    switch (sec) {
    case 1: da = t1 + t2 + half0; db = t2 + half0;  dc = half0;        break;
    case 2: da = t1 + half0;      db = t1 + t2 + half0; dc = half0;    break;
    case 3: da = half0;           db = t1 + t2 + half0; dc = t2 + half0; break;
    case 4: da = half0;           db = t1 + half0;  dc = t1 + t2 + half0; break;
    case 5: da = t2 + half0;      db = half0;       dc = t1 + t2 + half0; break;
    case 6: da = t1 + t2 + half0; db = half0;       dc = t1 + half0;   break;
    default: da = db = dc = 0.5; break;
    }
    out->ta = clamp01(da);
    out->tb = clamp01(db);
    out->tc = clamp01(dc);
}

/* ============================================================ 位置与速度
 * 霍尔换相表：120° 安装的三相霍尔，6 个有效状态。
 * 用 (hall-1)*60° 作为该区间的电角度基准（简化模型，实际需按安装角标定）。
 */
double sensor_hall_angle(uint8_t hall_raw, int *valid)
{
    const uint8_t h = hall_raw & 0x07u;
    /* 标准 120° 霍尔的 6 个有效状态 */
    static const uint8_t valid_states[6] = { 0x05, 0x01, 0x03, 0x02, 0x06, 0x04 };
    for (int i = 0; i < 6; i++) {
        if (h == valid_states[i]) {
            if (valid) *valid = 1;
            return (double)i * (FOC_PI / 3.0);
        }
    }
    if (valid) *valid = 0;            /* 000 与 111 是无效状态 */
    return 0.0;
}

void sensor_init(SensorState *s, double lpf_alpha)
{
    memset(s, 0, sizeof(*s));
    s->lpf_alpha = (lpf_alpha > 0.0 && lpf_alpha <= 1.0) ? lpf_alpha : 0.2;
}

void sensor_update_from_encoder(SensorState *s, int32_t count,
                                double dt, double pole_pairs)
{
    const int32_t delta = count - s->enc_last;
    s->enc_last = count;
    s->enc_count = count;

    /* 假设 1000 线编码器四倍频 -> 4000 counts/rev */
    const double cpr = 4000.0;
    s->theta_m = (double)count * (FOC_TWO_PI / cpr);
    s->theta_e = motor_elec_angle(s->theta_m, pole_pairs);

    if (dt > 0.0) {
        const double raw = (double)delta * (FOC_TWO_PI / cpr) / dt;
        s->omega_m = raw;
        /* 一阶低通，抑制差分噪声（低速时尤其需要） */
        s->omega_filt += s->lpf_alpha * (raw - s->omega_filt);
    }
    s->stamps++;
}

/* ========================================================= BLDC 六步换相 */
uint8_t bldc_next_step(uint8_t hall)
{
    static const uint8_t next[7] = { 0, 2, 3, 4, 5, 6, 1 };
    if (hall < 1 || hall > 6) return 0;
    return next[hall];
}

void bldc_six_step(uint8_t hall, double duty, SixStepOut *out)
{
    memset(out, 0, sizeof(*out));
    out->hall = hall;
    out->duty = duty < 0.0 ? 0.0 : (duty > 1.0 ? 1.0 : duty);
    out->step = hall;

    /* 标准六步换相：每步一相上桥、一相下桥、一相悬空（反电动势过零检测用） */
    static const uint8_t high[6] = { 0, 0, 1, 1, 2, 2 };
    static const uint8_t low[6]  = { 1, 2, 2, 0, 0, 1 };
    static const uint8_t floatp[6] = { 2, 1, 0, 2, 1, 0 };
    if (hall >= 1 && hall <= 6) {
        out->phase_high = high[hall - 1];
        out->phase_low  = low[hall - 1];
        out->floating   = floatp[hall - 1];
    }
}

/* ============================================================ 三相逆变器
 * 死区：上下桥臂切换时插入死区时间，防止直通短路。
 * 直通保护：同一相上下桥同时为 1 时置故障标志并计数。
 */
void inverter_init(InverterState *inv, double deadtime_us)
{
    memset(inv, 0, sizeof(*inv));
    inv->deadtime_us = deadtime_us;
}

void inverter_apply(InverterState *inv, const SvpwmOut *pwm, double udc)
{
    /* 上下桥臂**互补**：一个 PWM 周期内上管导通 ta 比例、下管导通 (1-ta)。
     * 任何瞬时只有一个导通，所以这里用互补表示，不会出现直通。
     * （要验证直通保护，可在单测里直接构造 bridge.bits 两位同时为 1） */
    inv->bridge.bits.ah = (pwm->ta > 0.5) ? 1 : 0;
    inv->bridge.bits.al = inv->bridge.bits.ah ? 0 : 1;
    inv->bridge.bits.bh = (pwm->tb > 0.5) ? 1 : 0;
    inv->bridge.bits.bl = inv->bridge.bits.bh ? 0 : 1;
    inv->bridge.bits.ch = (pwm->tc > 0.5) ? 1 : 0;
    inv->bridge.bits.cl = inv->bridge.bits.ch ? 0 : 1;
    inv->bridge.bits.deadtime_active = (inv->deadtime_us > 0.0) ? 1 : 0;

    /* 相电压（相对母线中点，简化模型） */
    inv->ua = (pwm->ta - 0.5) * udc;
    inv->ub = (pwm->tb - 0.5) * udc;
    inv->uc = (pwm->tc - 0.5) * udc;

    (void)inverter_check_shoot_through(inv);
}

int inverter_check_shoot_through(InverterState *inv)
{
    const uint16_t r = inv->bridge.raw;
    const int bad = ((r & 0x03u) == 0x03u) ||      /* A 相上下同通 */
                    ((r & 0x0Cu) == 0x0Cu) ||      /* B 相 */
                    ((r & 0x30u) == 0x30u);        /* C 相 */
    if (bad) {
        inv->bridge.bits.shoot_through = 1;
        inv->shoot_through_events++;
        return 1;
    }
    inv->bridge.bits.shoot_through = 0;
    return 0;
}

/* ================================================================== FOC */
FocConfig foc_default_config(void)
{
    FocConfig c;
    c.kp_id = 8.0;   c.ki_id = 600.0;
    c.kp_iq = 8.0;   c.ki_iq = 600.0;
    c.kp_sp = 0.06;  c.ki_sp = 0.6;
    c.ts = 1.0 / 20000.0;      /* 20 kHz 电流环，车用控制器常见 */
    c.udc = 48.0;      /* 母线 48 V：24 V 下反电动势会先顶到电压极限 */
    c.i_max = 20.0;
    c.u_max_ratio = 0.95;
    c.speed_lpf_alpha = 0.25;
    return c;
}

void foc_init(FocCtrl *c, const FocConfig *cfg)
{
    memset(c, 0, sizeof(*c));
    const double umax = cfg->udc * FOC_1_SQRT3 * cfg->u_max_ratio;
    pi_init(&c->id_pi, cfg->kp_id, cfg->ki_id, cfg->ts, -umax, umax);
    pi_init(&c->iq_pi, cfg->kp_iq, cfg->ki_iq, cfg->ts, -umax, umax);
    pi_init(&c->speed_pi, cfg->kp_sp, cfg->ki_sp, cfg->ts,
            -cfg->i_max, cfg->i_max);
    c->udc = cfg->udc;
    c->i_ref.d = 0.0;          /* id = 0 控制（表贴式 PMSM 的最优转矩/电流比） */
    c->i_ref.q = 0.0;
}

double foc_field_weakening(double omega_e, double omega_base, double id_max)
{
    if (omega_base <= 0.0 || omega_e <= omega_base) return 0.0;
    /* 转速超过基速后按比例给负 id，削弱磁链以扩速 */
    const double ratio = (omega_e - omega_base) / omega_base;
    double id = -id_max * ratio;
    if (id < -id_max) id = -id_max;
    return id;
}

void foc_step(FocCtrl *c, const FocConfig *cfg, const MotorParams *p,
              MotorState *motor, SensorState *sensor, InverterState *inv,
              double dt)
{
    c->speed_fb = sensor->omega_filt;

    /* 1) 速度环 -> iq 参考（带转矩前馈） */
    pi_set_feedforward(&c->speed_pi,
                       c->torque_ff / (1.5 * p->pole_pairs * p->flux + 1e-12));
    c->i_ref.q = pi_step(&c->speed_pi, c->speed_ref, c->speed_fb);

    /* 2) 电流环（id / iq 各自 PI，反馈来自电机模型） */
    c->i_fb.d = motor->id;
    c->i_fb.q = motor->iq;
    c->u_ref.d = pi_step(&c->id_pi, c->i_ref.d, c->i_fb.d);
    c->u_ref.q = pi_step(&c->iq_pi, c->i_ref.q, c->i_fb.q);

    /* 3) 前馈解耦：补偿 dq 轴的交叉耦合与反电动势
     *    ud_ff = -ωe·Lq·iq , uq_ff = ωe·(Ld·id + ψf)
     *    加上前馈后电流环在高转速下仍能保持带宽。 */
    const double we = motor->omega_e;
    c->u_ref.d += -we * p->lq * c->i_fb.q;
    c->u_ref.q +=  we * (p->ld * c->i_fb.d + p->flux);

    /* 4) 电压限幅 */
    const double umax = cfg->udc * FOC_1_SQRT3 * cfg->u_max_ratio;
    const double umag = sqrt(c->u_ref.d * c->u_ref.d + c->u_ref.q * c->u_ref.q);
    if (umag > umax && umag > 0.0) {
        const double s = umax / umag;
        c->u_ref.d *= s;
        c->u_ref.q *= s;
    }

    /* 5) SVPWM -> 三相占空比 -> 驱动逆变器 -> 作用到电机 */
    svpwm_calc(c->u_ref, motor->theta_e, cfg->udc, &c->pwm);
    inverter_apply(inv, &c->pwm, cfg->udc);

    /* 把相电压变换回 dq 施加到电机模型（说明：真实系统由桥臂电压产生） */
    const Abc uabc = { inv->ua, inv->ub, inv->uc };
    const Dq udq = park(clarke(uabc), motor->theta_e);
    motor->ud = udq.d;
    motor->uq = udq.q;

    c->isr_count++;
}

/* ================================================================== 统计 */
void foc_stats_reset(FocStats *s)
{
    memset(s, 0, sizeof(*s));
}

void foc_stats_update(FocStats *s, const MotorState *m, const MotorParams *p,
                      const InverterState *inv, double udc)
{
    const double ia = m->ia, ib = m->ib, ic = m->ic;
    /* 增量式均值（避免长期累加溢出） */
    const double n = (double)(s->samples + 1);
    s->ia_rms += (ia * ia - s->ia_rms) / n;
    s->ib_rms += (ib * ib - s->ib_rms) / n;
    s->ic_rms += (ic * ic - s->ic_rms) / n;
    s->id_mean += (m->id - s->id_mean) / n;
    s->iq_mean += (m->iq - s->iq_mean) / n;
    s->speed_mean += (m->omega_m - s->speed_mean) / n;
    s->speed_ripple += (fabs(m->omega_m - s->speed_mean) - s->speed_ripple) / n;
    s->torque_mean += (m->te - s->torque_mean) / n;
    s->torque_ripple += (fabs(m->te - s->torque_mean) - s->torque_ripple) / n;

    /* 输入功率用 dq 轴口径：P = 1.5·(ud·id + uq·iq)，而不是相电流绝对值之和 */
    const double p_in = 1.5 * (m->ud * m->id + m->uq * m->iq);
    const double p_out = m->te * m->omega_m;
    if (fabs(p_in) > 1e-6) {
        /* 只在电动状态（输入为正且输出为正）统计效率，避免回馈/空载拉低均值 */
        if (p_in > 1e-3 && p_out > 0.0) {
            const double eff = p_out / p_in;
            s->efficiency += (eff - s->efficiency) / n;
        }
    }
    (void)ia; (void)ib; (void)ic; (void)udc;
    s->shoot_through = inv ? inv->shoot_through_events : 0;
    s->samples++;
    (void)p;
}
