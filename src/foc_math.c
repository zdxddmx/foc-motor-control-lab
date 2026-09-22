/* foc_math.c — 电机模型、坐标变换、PI 控制器 */
#include "foc.h"
#include <math.h>
#include <string.h>

MotorParams motor_default_pmsm(void)
{
    MotorParams p;
    p.rs = 0.35;            /* Ω    */
    p.ld = 0.0008;          /* H    */
    p.lq = 0.0008;          /* H（表贴式：Ld=Lq） */
    p.flux = 0.0175;        /* Wb   */
    p.pole_pairs = 4.0;
    p.j = 1.2e-5;           /* kg·m² */
    p.b = 2.0e-6;
    p.ke = 0.0175 * 1.5 * 4.0;
    return p;
}

MotorParams motor_default_bldc(void)
{
    MotorParams p = motor_default_pmsm();
    p.ld = 0.0006;          /* 无刷直流常见凸极：Ld != Lq */
    p.lq = 0.0009;
    p.flux = 0.0200;
    p.ke = 0.0200 * 1.5 * 4.0;
    return p;
}

void motor_init(MotorState *st)
{
    memset(st, 0, sizeof(*st));
}

double motor_elec_angle(double theta_m, double pole_pairs)
{
    double e = theta_m * pole_pairs;
    e = fmod(e, FOC_TWO_PI);
    if (e < 0.0) e += FOC_TWO_PI;
    return e;
}

/* ------------------------------------------------------------ 电机模型
 * dq 轴方程（表贴式 PMSM，Ld=Lq，含反电动势）：
 *   did/dt = (ud - Rs·id + ωe·Lq·iq) / Ld
 *   diq/dt = (uq - Rs·iq - ωe·Ld·id - ωe·ψf) / Lq
 *   电磁转矩 Te = 1.5·P·[ψf·iq + (Ld-Lq)·id·iq]
 *   机械方程 J·dωm/dt = Te - TL - B·ωm
 */
void motor_step(MotorState *st, const MotorParams *p, double dt)
{
    const double we = st->omega_e;
    const double did = (st->ud - p->rs * st->id + we * p->lq * st->iq) / p->ld;
    const double diq = (st->uq - p->rs * st->iq - we * p->ld * st->id - we * p->flux)
                       / p->lq;

    st->id += did * dt;
    st->iq += diq * dt;

    st->te = 1.5 * p->pole_pairs *
             (p->flux * st->iq + (p->ld - p->lq) * st->id * st->iq);

    const double dw = (st->te - st->tl - p->b * st->omega_m) / p->j;
    st->omega_m += dw * dt;
    st->omega_e = st->omega_m * p->pole_pairs;
    st->theta_e = motor_elec_angle(st->theta_e + st->omega_e * dt, 1.0);

    /* 反推三相电流（用于验证变换往返与实际相电流） */
    const AlphaBeta ab = inv_park((Dq){ st->id, st->iq }, st->theta_e);
    const Abc abc = inv_clarke(ab);
    st->ia = abc.a;
    st->ib = abc.b;
    st->ic = abc.c;
    st->steps++;
}

void motor_abc_to_dq(const double iabc[3], double theta_e, double *id, double *iq)
{
    const Abc in = { iabc[0], iabc[1], iabc[2] };
    const Dq out = park(clarke(in), theta_e);
    *id = out.d;
    *iq = out.q;
}

/* ---------------------------------------------------------------- 变换
 * 等幅值 Clarke（与常见工程实现一致）：
 *   α = (2/3)(ia - ib/2 - ic/2)
 *   β = (2/3)(√3/2)(ib - ic) = (ib - ic)/√3
 */
AlphaBeta clarke(Abc in)
{
    AlphaBeta o;
    o.alpha = (2.0 / 3.0) * (in.a - 0.5 * in.b - 0.5 * in.c);
    o.beta = (in.b - in.c) * FOC_1_SQRT3;
    return o;
}

Abc inv_clarke(AlphaBeta in)
{
    Abc o;
    o.a = in.alpha;
    o.b = -0.5 * in.alpha + FOC_SQRT3_2 * in.beta;
    o.c = -0.5 * in.alpha - FOC_SQRT3_2 * in.beta;
    return o;
}

Dq park(AlphaBeta in, double theta_e)
{
    const double c = cos(theta_e), s = sin(theta_e);
    Dq o;
    o.d = in.alpha * c + in.beta * s;
    o.q = -in.alpha * s + in.beta * c;
    return o;
}

AlphaBeta inv_park(Dq in, double theta_e)
{
    const double c = cos(theta_e), s = sin(theta_e);
    AlphaBeta o;
    o.alpha = in.d * c - in.q * s;
    o.beta = in.d * s + in.q * c;
    return o;
}

/* ---------------------------------------------------------------- PI
 * 带条件积分的抗饱和：
 *   只有"输出未饱和"或"误差方向会把输出拉回线性区"时才累加积分项，
 *   避免积分饱和导致的超调与恢复迟滞（电机启动/堵转时很关键）。
 */
void pi_init(PiCtrl *pi, double kp, double ki, double ts,
             double out_min, double out_max)
{
    memset(pi, 0, sizeof(*pi));
    pi->kp = kp;
    pi->ki = ki;
    pi->ts = ts;
    pi->out_min = out_min;
    pi->out_max = out_max;
}

void pi_set_feedforward(PiCtrl *pi, double ff)
{
    pi->ff = ff;
}

double pi_step(PiCtrl *pi, double ref, double fb)
{
    const double err = ref - fb;
    double p = pi->kp * err;
    double i_term = pi->ki * pi->integ;
    double out = p + i_term + pi->ff;

    /* 条件积分：饱和且误差继续把输出推向饱和方向时，不累加 */
    const int push_high = (out > pi->out_max) && (err > 0.0);
    const int push_low = (out < pi->out_min) && (err < 0.0);
    if (!push_high && !push_low) {
        pi->integ += err * pi->ts;
        i_term = pi->ki * pi->integ;
        out = p + i_term + pi->ff;
    }
    pi->saturated = (out >= pi->out_max) || (out <= pi->out_min);

    if (out > pi->out_max) out = pi->out_max;
    if (out < pi->out_min) out = pi->out_min;
    pi->last_err = err;
    pi->calls++;
    return out;
}

void pi_reset(PiCtrl *pi)
{
    pi->integ = 0.0;
    pi->last_err = 0.0;
    pi->saturated = 0;
}
