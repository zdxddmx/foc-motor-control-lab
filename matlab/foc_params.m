function p = foc_params()
% FOC_PARAMS  电机与控制器参数（含 PI 整定）
%
% 用途：把"参数从哪来"讲清楚，并给出一套可核对的整定结果。
%   1) 电机电气时间常数 tau = L/R  -> 决定电流环带宽上限
%   2) 电流环按"期望带宽"整定：Kp = L * wc, Ki = R * wc
%      （把被控对象近似为一阶 L/R 环节，PI 零点抵消极点）
%   3) 速度环带宽取电流环的 1/10（带宽分离原则），
%      Kp = J * wv / (1.5 * P * psi_f), Ki 取 Kp * wv / 10
%   4) 母线电压 -> 可达到的最大电角速度（弱磁起点）
%
% 兼容 GNU Octave 与 MATLAB，不依赖任何工具箱。
%
% 用法：
%   p = foc_params(); disp(p)
%   p = foc_params(); p.plot = true; foc_params(p);   % 可选画图

  if nargin < 1
    plot_figs = false;
  else
    plot_figs = false;
  end

  %% ---------- 电机参数（表贴式 PMSM，与 C 版 motor_default_pmsm 一致）----------
  p.Rs        = 0.35;        % 定子电阻 [Ohm]
  p.Ld        = 0.8e-3;      % d 轴电感 [H]
  p.Lq        = 0.8e-3;      % q 轴电感 [H]（表贴式：Ld = Lq）
  p.psi_f     = 0.0175;      % 永磁磁链 [Wb]
  p.pole_pairs= 4;           % 极对数
  p.J         = 1.2e-5;      % 转动惯量 [kg*m^2]
  p.B         = 2.0e-6;      % 粘滞摩擦 [N*m*s/rad]

  %% ---------- 逆变器与采样 ----------
  p.Udc       = 48;          % 母线电压 [V]
  p.f_sw      = 20000;       % 开关/电流环频率 [Hz]
  p.Ts        = 1 / p.f_sw;  % 电流环周期 [s]
  p.deadtime  = 1e-6;        % 死区 [s]

  %% ---------- 电气时间常数与电流环整定 ----------
  p.tau_e     = p.Ld / p.Rs;                 % 电气时间常数 [s]
  p.wc_target = 2 * pi * 800;                % 期望电流环带宽 [rad/s]（800 Hz）
  % 一阶对象 L/R，PI 零点抵消极点
  p.Kp_i      = p.Ld * p.wc_target;
  p.Ki_i      = p.Rs * p.wc_target;
  % 校验：带宽不能超过开关频率的 1/10（否则离散化误差大）
  p.wc_limit  = 2 * pi * p.f_sw / 10;
  p.wc_ok     = p.wc_target <= p.wc_limit;

  %% ---------- 电压极限与弱磁起点 ----------
  p.Umax      = p.Udc / sqrt(3) * 0.95;      % SVPWM 线性区电压极限 [V]
  p.we_base   = p.Umax / p.psi_f;            % 基速（电角速度）[rad/s]
  p.wm_base   = p.we_base / p.pole_pairs;    % 基速（机械角速度）[rad/s]
  p.imax      = 20;                          % 电流限幅 [A]

  %% ---------- 速度环整定（带宽分离：取电流环 1/10）----------
  p.wv_target = p.wc_target / 10;            % 速度环带宽 [rad/s]
  p.Kt        = 1.5 * p.pole_pairs * p.psi_f;% 转矩常数 [N*m/A]
  p.Kp_w      = p.J * p.wv_target / p.Kt;
  p.Ki_w      = p.Kp_w * p.wv_target / 10;

  %% ---------- 打印 ----------
  fprintf('==== FOC 参数与整定 ====\n');
  fprintf('电机：Rs=%.3f Ohm  Ld=Lq=%.4f mH  psi_f=%.4f Wb  P=%d\n', ...
          p.Rs, p.Ld * 1e3, p.psi_f, p.pole_pairs);
  fprintf('电气时间常数 tau_e = %.3f ms\n', p.tau_e * 1e3);
  fprintf('电流环：目标带宽 %.0f Hz，Kp=%.5f  Ki=%.3f（上限 %.0f Hz，%s）\n', ...
          p.wc_target / (2*pi), p.Kp_i, p.Ki_i, p.wc_limit / (2*pi), ...
          merge(p.wc_ok, '通过', '超出'));
  fprintf('速度环：目标带宽 %.0f Hz，Kp=%.5f  Ki=%.5f  转矩常数 Kt=%.5f N*m/A\n', ...
          p.wv_target / (2*pi), p.Kp_w, p.Ki_w, p.Kt);
  fprintf('电压极限 Umax = %.2f V  ->  基速 %.0f rad/s（电）/ %.0f rad/s（机械）\n', ...
          p.Umax, p.we_base, p.wm_base);

  %% ---------- 可选画图（需图形环境；CI 下默认关闭）----------
  if plot_figs
    w = logspace(1, 5, 400);
    % 电流环开环：PI * 一阶对象
    s = 1i * w;
    G = (p.Kp_i + p.Ki_i ./ s) .* (1 ./ (p.Rs + p.Ld .* s));
    figure; semilogx(w / (2*pi), 20*log10(abs(G))); grid on;
    xlabel('Hz'); ylabel('|G| dB'); title('电流环开环幅频');
    print('-dpng', 'foc_current_loop_bode.png');
  end
end

function out = merge(cond, a, b)
  if cond, out = a; else, out = b; end
end
