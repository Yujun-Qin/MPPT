% PV 单二极管模型 + 理想 Boost 平均模型 + 变步长(基于 dP/dV) P&O MPPT
% 运行环境：MATLAB/Octave (fzero在Octave需optim包)

% ==== PV 组件与物理参数（可按实际组件微调）====
pv.Ns        = 60;      % 串联片数
pv.Rs        = 0.35;    % 串联电阻 [Ohm]
pv.Rsh       = 500;     % 并联电阻 [Ohm]
pv.n         = 1.3;     % 理想因子
pv.Isc_stc   = 8.7;     % STC短路电流 [A]
pv.Voc_stc   = 37.0;    % STC开路电压 [V]
pv.alpha_Isc = 0.005;   % A/°C
pv.Eg        = 1.12;    % eV

G  = 800;               % W/m^2  (环境辐照)
Tc = 35;                % °C     (组件温度)

% ==== 负载与Boost平均模型 ====
Rload = 20;             % 直流负载 [Ohm]（等效纯电阻）
Dmin = 0.02; Dmax = 0.95;   % 占空比边界（避免极端）
D    = 0.5;                 % 初始占空比

% 说明：在理想连续导通(CCM)平均模型下，源侧看到的等效阻抗
% Rin = Rload * (1 - D)^2
% 工作点满足 Ipv(V) = V / Rin。给定D可解得 Vpv。

% ==== 展示静态 I–V / P–V（当前G,T）以及MPPT最终点 ====
[Vvec, Ivec] = sweep_pv_iv(G, Tc, pv);
Pvec = Vvec .* Ivec;
V_max = max(Vvec);
I_max = max()
P_max = max(Pvec);

figure;
plot(Vvec, Ivec, 'LineWidth',1.5); grid on; hold on;
xlabel('V_{pv} (V)'); ylabel('I_{pv} (A)');
title('PV I_V (current G,T) & MPPT result');
legend('I_V','MPPT point','Location','best');

figure;
plot(Vvec, Pvec, 'LineWidth',1.5); grid on; hold on;

xlabel('V_{pv} (V)'); ylabel('P_{pv} (W)');
title('PV P_V (current G,T) & MPPT result');

% ===== 局部函数 =====
function [Vsol, Isol] = solve_pv_operating_point(D, Rload, G, T, params)
    % 给定占空比D与负载Rload，在理想Boost平均模型下解方程：
    % Ipv(V) = V / Rin，其中 Rin = Rload*(1-D)^2
    Rin = Rload * (1 - D)^2;
    if Rin <= 0, Rin = 1e-6; end

    % 残差方程： f(V) = Ipv(V) - V/Rin = 0
    f = @(V) pv_current_given_voltage(max(V,0), G, T, params) - max(V,0)/Rin;

    % 解V：初值与区间（在[0, Voc估计]内）
    Voc_est = params.Voc_stc * (1 + 0.0*(T-25)/25); % 简化：Voc弱温度依赖
    V0 = 0.7*Voc_est;
    try
        Vsol = fzero(f, [0, max(1e-6, 1.2*Voc_est)]);
    catch
        Vsol = fzero(f, V0);
    end
    Vsol = max(Vsol, 0);
    Isol = pv_current_given_voltage(Vsol, G, T, params);
end

function I = pv_current_given_voltage(V, G, T, params)
    % 单二极管模型数值解 I(V):
    % I = Iph - I0*(exp((V+I*Rs)/(n*Vt)) - 1) - (V+I*Rs)/Rsh
    q = 1.602176634e-19; k_b = 1.380649e-23;
    T_K = T + 273.15; 
    Vt  = params.Ns * k_b*T_K / q;
    n   = params.n;

    % 光生电流 Iph
    Iph_stc = params.Isc_stc;
    Iph = Iph_stc*(G/1000) + params.alpha_Isc*(T - 25);

    % 反向饱和电流 I0 温度修正（由Voc_stc反推I0_stc再修正）
    Vt_stc = params.Ns * k_b*(25+273.15)/q;
    I0_stc = Iph_stc / (exp(params.Voc_stc/(n*Vt_stc)) - 1 + 1e-12);
    Eg = params.Eg;
    I0 = I0_stc * (T_K/(25+273.15)).^3 .* ...
         exp( (q*Eg/k_b) * (1/(25+273.15) - 1/T_K) / (n*params.Ns) );

    Rs  = params.Rs; Rsh = params.Rsh;

    fI = @(I) Iph - I0.*(exp((V + I.*Rs)./(n*Vt)) - 1) - (V + I.*Rs)./Rsh - I;
    I_init = max(0, min(Iph, Iph - V/max(Rsh,1e-6)));
    try
        I = fzero(fI, I_init);
    catch
        I = fzero(fI, [0, max(0.1, Iph*1.5)]);
    end
    I = max(I, 0);
end

function [Vvec, Ivec] = sweep_pv_iv(G, T, params)
    % 生成当前G,T下的I–V曲线，便于可视化
    Vmax = params.Voc_stc*1.15;
    Vvec = linspace(0, Vmax, 200);
    Ivec = zeros(size(Vvec));
    for i = 1:numel(Vvec)
        Ivec(i) = pv_current_given_voltage(Vvec(i), G, T, params);
    end
end
