function out = foc_svpwm(ualpha, ubeta, Udc)
% FOC_SVPWM  空间矢量脉宽调制（与 C 版 svpwm_calc 同一套算法）
%
% 输入：ualpha/ubeta（静止坐标系电压，V），Udc（母线电压，V）
% 输出：结构体，含扇区、T1/T2/T0、三相占空比、调制比
%
% 用于核对 C 实现的正确性：同样的输入应给出同样的占空比。
% 兼容 GNU Octave 与 MATLAB，无工具箱依赖。

  if nargin < 3
    Udc = 48;
  end

  out.ualpha = ualpha;
  out.ubeta  = ubeta;

  umag = sqrt(ualpha^2 + ubeta^2);
  if Udc > 0
    out.mod_index = umag / (Udc / sqrt(3));
  else
    out.mod_index = 0;
  end

  %% ---- 扇区判定：用三个判断量，避免 atan2（车用 MCU 上更省时间）----
  a = ubeta;
  b = -0.5 * ubeta + sqrt(3)/2 * ualpha;
  c = -0.5 * ubeta - sqrt(3)/2 * ualpha;
  s = 0;
  if a > 0, s = s + 1; end
  if b > 0, s = s + 2; end
  if c > 0, s = s + 4; end
  switch s
    case 3, sec = 1;
    case 1, sec = 2;
    case 5, sec = 3;
    case 4, sec = 4;
    case 6, sec = 5;
    case 2, sec = 6;
    otherwise, sec = 0;      % 零矢量附近
  end
  out.sector = sec;

  if sec == 0 || Udc <= 0
    out.T1 = 0; out.T2 = 0; out.T0 = 1;
    out.ta = 0.5; out.tb = 0.5; out.tc = 0.5;
    return;
  end

  %% ---- 相邻基本矢量作用时间 ----
  k = sqrt(3) / Udc;
  switch sec
    case 1
      T1 = k * ( sqrt(3)/2 * ualpha - 0.5 * ubeta);
      T2 = k * ( sqrt(3)/2 * ualpha + 0.5 * ubeta);
    case 2
      T1 = k * ( sqrt(3)/2 * ualpha + 0.5 * ubeta);
      T2 = k * (-sqrt(3)/2 * ualpha + 0.5 * ubeta);
    case 3
      T1 = k * ( ubeta);
      T2 = k * (-ualpha);
    case 4
      T1 = k * (-sqrt(3)/2 * ualpha + 0.5 * ubeta);
      T2 = k * (-sqrt(3)/2 * ualpha - 0.5 * ubeta);
    case 5
      T1 = k * (-sqrt(3)/2 * ualpha - 0.5 * ubeta);
      T2 = k * ( sqrt(3)/2 * ualpha + 0.5 * ubeta);
    case 6
      T1 = k * (-ubeta);
      T2 = k * ( ualpha);
    otherwise
      T1 = 0; T2 = 0;
  end
  if T1 < 0, T1 = 0; end
  if T2 < 0, T2 = 0; end

  %% ---- 过调制保护：T1 + T2 不得超过 1 ----
  if T1 + T2 > 1
    sc = 1 / (T1 + T2);
    T1 = T1 * sc;
    T2 = T2 * sc;
    out.overmodulated = true;
  else
    out.overmodulated = false;
  end
  T0 = 1 - T1 - T2;
  out.T1 = T1; out.T2 = T2; out.T0 = T0;

  %% ---- 五段式：把 T1/T2 按扇区分配到三相 ----
  h = T0 * 0.5;
  switch sec
    case 1, da = T1+T2+h; db = T2+h;    dc = h;
    case 2, da = T1+h;    db = T1+T2+h; dc = h;
    case 3, da = h;       db = T1+T2+h; dc = T2+h;
    case 4, da = h;       db = T1+h;    dc = T1+T2+h;
    case 5, da = T2+h;    db = h;       dc = T1+T2+h;
    case 6, da = T1+T2+h; db = h;       dc = T1+h;
    otherwise, da = 0.5; db = 0.5; dc = 0.5;
  end
  out.ta = min(max(da, 0), 1);
  out.tb = min(max(db, 0), 1);
  out.tc = min(max(dc, 0), 1);

  %% ---- 线电压校验：三者之和应为 0 ----
  out.uab = (out.ta - out.tb) * Udc;
  out.ubc = (out.tb - out.tc) * Udc;
  out.uca = (out.tc - out.ta) * Udc;
end
