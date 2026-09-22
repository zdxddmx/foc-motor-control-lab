/*
 * foc.h — 无刷电机 FOC 控制库（BLDC / PMSM）
 *
 * 对齐汽车电子电机控制器岗位的技术栈：
 *   三相逆变器 + SVPWM       -> svpwm.c
 *   Clarke / Park 变换       -> foc.c
 *   电流环 / 速度环 PI + 前馈 -> foc.c
 *   BLDC 六步换相 vs PMSM FOC -> bldc.c
 *   霍尔 / 编码器位置与速度   -> sensor.c
 *   电机模型（dq 方程 + 反电动势）-> motor.c
 *
 * 代码风格刻意贴合车用 MCU 固件的写法（JD 要求"结构体、共用体、指针、内存占用"）：
 *   · 定点/浮点混用的数据结构用 struct 收口，避免全局散变量
 *   · 寄存器位域与状态用 union 表达（同一块内存两种视图）
 *   · 运行时数据用指针传入，函数不依赖全局状态，便于上板与单测
 *   · 明确标注内存占用，便于评估 RAM 预算
 *
 * 纯 C99，零第三方依赖。
 */
#ifndef FOC_H
#define FOC_H

#include <stdint.h>
#include <stddef.h>

/* ---------------------------------------------------------------- 数学常量 */
#define FOC_PI        3.14159265358979323846
#define FOC_TWO_PI    (2.0 * FOC_PI)
#define FOC_SQRT3     1.73205080756887729352
#define FOC_1_SQRT3   0.57735026918962576451
#define FOC_SQRT3_2   0.86602540378443864676

/* ---------------------------------------------------------------- 电机参数 */
typedef struct {
    double rs;          /* 定子电阻 Ω */
    double ld;          /* d 轴电感 H */
    double lq;          /* q 轴电感 H */
    double flux;        /* 永磁磁链 Wb */
    double pole_pairs;  /* 极对数 */
    double j;           /* 转动惯量 kg·m² */
    double b;           /* 粘滞摩擦系数 */
    double ke;          /* 反电动势常数 V·s/rad（表贴式时 = flux * 1.5 * P） */
} MotorParams;

/* 表贴式 PMSM：Ld = Lq */
MotorParams motor_default_pmsm(void);
MotorParams motor_default_bldc(void);

/* ------------------------------------------------------------ 电机状态/模型 */
typedef struct {
    double id, iq;              /* dq 轴电流 A */
    double ia, ib, ic;          /* 三相电流 A（用于验证变换往返） */
    double theta_e;             /* 电角度 rad */
    double omega_m;             /* 机械角速度 rad/s */
    double omega_e;             /* 电角速度 rad/s */
    double te;                  /* 电磁转矩 N·m */
    double tl;                  /* 负载转矩 N·m */
    double ud, uq;              /* 施加的 dq 电压 V */
    uint64_t steps;             /* 仿真步数 */
} MotorState;

void motor_init(MotorState *st);
/* 按 dq 方程推进一个步长：电机模型（含反电动势与转矩方程） */
void motor_step(MotorState *st, const MotorParams *p, double dt);
/* 机械角速度 -> 电角度（并做 0~2π 归一化） */
double motor_elec_angle(double theta_m, double pole_pairs);
/* 三相电流 -> dq（Clarke + Park），用于校验 */
void motor_abc_to_dq(const double iabc[3], double theta_e, double *id, double *iq);

/* ---------------------------------------------------------------- 变换 */
typedef struct { double a, b, c; } Abc;
typedef struct { double alpha, beta; } AlphaBeta;
typedef struct { double d, q; } Dq;

AlphaBeta clarke(Abc in);                 /* 等幅值 Clarke */
Abc       inv_clarke(AlphaBeta in);
Dq        park(AlphaBeta in, double theta_e);
AlphaBeta inv_park(Dq in, double theta_e);

/* ---------------------------------------------------------------- PI 控制器 */
typedef struct {
    double kp, ki;
    double ts;              /* 采样周期 s */
    double out_min, out_max;
    double integ;           /* 积分累加（带抗饱和） */
    double ff;              /* 前馈项（JD 要求"前馈控制"） */
    double last_err;
    int    saturated;       /* 上次是否饱和（条件积分用） */
    uint32_t calls;
} PiCtrl;

void   pi_init(PiCtrl *pi, double kp, double ki, double ts,
               double out_min, double out_max);
void   pi_set_feedforward(PiCtrl *pi, double ff);
double pi_step(PiCtrl *pi, double ref, double fb);
void   pi_reset(PiCtrl *pi);

/* ------------------------------------------------------------------ SVPWM */
typedef struct {
    double ta, tb, tc;      /* 三相占空比 0~1 */
    double t1, t2;          /* 两个基本矢量作用时间（归一化到 Ts） */
    double t0;              /* 零矢量时间 */
    uint8_t sector;         /* 扇区 1~6，0 表示零矢量 */
    double mod_index;       /* 调制比 */
} SvpwmOut;

/* 输入 dq 电压与电角度，输出三相占空比；母线电压 udc 用于限幅 */
void svpwm_calc(Dq udq, double theta_e, double udc, SvpwmOut *out);
/* 扇区判定（独立出来便于单测） */
uint8_t svpwm_sector(double alpha, double beta);

/* ------------------------------------------------------ 位置/速度传感器 */
typedef struct {
    /* union：同一块内存两种视图 —— 车用固件里常见的写法 */
    union {
        struct {
            uint8_t h1 : 1;
            uint8_t h2 : 1;
            uint8_t h3 : 1;
            uint8_t dir : 1;      /* 转向 */
            uint8_t fault : 1;
            uint8_t : 3;          /* 补齐 1 字节 */
        } bits;
        uint8_t raw;
    } hall;
    int32_t  enc_count;      /* 编码器计数 */
    int32_t  enc_last;
    double   theta_m;        /* 机械角 rad */
    double   theta_e;        /* 电角 rad */
    double   omega_m;        /* 机械角速度 rad/s（差分 + 低通） */
    double   omega_filt;
    double   lpf_alpha;      /* 速度低通系数 */
    uint64_t stamps;
} SensorState;

void   sensor_init(SensorState *s, double lpf_alpha);
/* 用霍尔状态估算电角度（六步换相用）；返回该霍尔区间对应的角度 */
double sensor_hall_angle(uint8_t hall_raw, int *valid);
/* 由编码器计数更新角度与速度 */
void   sensor_update_from_encoder(SensorState *s, int32_t count,
                                  double dt, double pole_pairs);

/* ------------------------------------------------------- BLDC 六步换相 */
typedef struct {
    uint8_t step;           /* 1~6 */
    uint8_t hall;           /* 当前霍尔状态 */
    double  duty;           /* 当前占空比 */
    uint8_t phase_high;     /* 上桥臂导通的相 0=A,1=B,2=C */
    uint8_t phase_low;      /* 下桥臂导通的相 */
    uint8_t floating;       /* 悬空相（用于反电动势过零检测） */
} SixStepOut;

/* 标准 120° 霍尔换相表：hall(1~6) -> 下一步 */
uint8_t bldc_next_step(uint8_t hall);
void    bldc_six_step(uint8_t hall, double duty, SixStepOut *out);

/* ------------------------------------------------------ 三相逆变器桥 */
typedef struct {
    /* 上下桥臂开关状态：同一相不允许直通（死区由 deadtime 处理） */
    union {
        struct {
            uint16_t ah : 1, al : 1;
            uint16_t bh : 1, bl : 1;
            uint16_t ch : 1, cl : 1;
            uint16_t deadtime_active : 1;
            uint16_t shoot_through : 1;   /* 直通故障标志 */
            uint16_t : 8;
        } bits;
        uint16_t raw;
    } bridge;
    double deadtime_us;      /* 死区时间 us */
    double ua, ub, uc;       /* 相电压 V */
    uint32_t shoot_through_events;
} InverterState;

void inverter_init(InverterState *inv, double deadtime_us);
/* 由三相占空比驱动桥臂（带死区与直通保护） */
void inverter_apply(InverterState *inv, const SvpwmOut *pwm, double udc);
/* 检查是否发生上下桥直通 */
int  inverter_check_shoot_through(InverterState *inv);

/* ------------------------------------------------------------------ 控制器 */
typedef struct {
    PiCtrl  id_pi;          /* d 轴电流环 */
    PiCtrl  iq_pi;          /* q 轴电流环 */
    PiCtrl  speed_pi;       /* 速度环 */
    Dq      i_ref;          /* 电流参考（id 通常置 0，即 id=0 控制） */
    Dq      i_fb;           /* 电流反馈 */
    Dq      u_ref;          /* 电压参考 */
    double  speed_ref;      /* 机械角速度参考 rad/s */
    double  speed_fb;
    double  torque_ff;      /* 转矩前馈（负载已知时直接前馈） */
    SvpwmOut pwm;
    double  udc;
    uint64_t isr_count;
} FocCtrl;

typedef struct {
    double kp_id, ki_id;
    double kp_iq, ki_iq;
    double kp_sp, ki_sp;
    double ts;              /* 电流环周期 s */
    double udc;
    double i_max;           /* 电流限幅 A */
    double u_max_ratio;     /* 电压限幅（占母线比例） */
    double speed_lpf_alpha;
} FocConfig;

FocConfig foc_default_config(void);
void foc_init(FocCtrl *c, const FocConfig *cfg);
/* 一个电流环周期：采样 -> PI -> 解耦/前馈 -> SVPWM */
void foc_step(FocCtrl *c, const FocConfig *cfg, const MotorParams *p,
              MotorState *motor, SensorState *sensor, InverterState *inv,
              double dt);
/* 弱磁（高速时给负 id）：可选 */
double foc_field_weakening(double omega_e, double omega_base, double id_max);

/* ------------------------------------------------------------------ 统计 */
typedef struct {
    double   ia_rms, ib_rms, ic_rms;
    double   id_mean, iq_mean;
    double   speed_mean, speed_ripple;
    double   torque_mean, torque_ripple;
    double   efficiency;        /* 机械输出 / 电输入 */
    uint64_t samples;
    uint32_t shoot_through;
} FocStats;

void foc_stats_reset(FocStats *s);
void foc_stats_update(FocStats *s, const MotorState *m, const MotorParams *p,
                      const InverterState *inv, double udc);

#endif /* FOC_H */
