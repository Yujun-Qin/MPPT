function [Plog] = mppt_INC()
% pv_mppt_boost_slopePO.m
% 单文件：PV 单二极管模型 + 理想 Boost 平均模型 + 变步长(基于 dP/dV) P&O MPPT
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

% ==== MPPT：增量电导法（Incremental Conductance，基于探测点的实现）====

Nstep    = 100;      % 迭代步数
DeltaD   = 0.001;     % 占空比固定步长（决定靠近 MPP 的速度&抖动）
test_step = 1e-4;    % 用于估计 dI/dV 的极小扰动 (在D维度里)

epsV    = 1e-4;      % 电压判零阈值
epsI    = 1e-4;      % 电流判零阈值
epsCond = 1e-3;      % dI/dV + I/V 接近0的判定阈值

% 记录量
Dlog = zeros(1,Nstep);
Vlog = zeros(1,Nstep);
Ilog = zeros(1,Nstep);
Plog = zeros(1,Nstep);

% 初始工作点 (使用初始占空比 D)
[Vk, Ik] = solve_pv_operating_point(D, Rload, G, Tc, pv);

for k = 1:Nstep
    % 1) 记录当前工作点
    Dlog(k) = D;
    Vlog(k) = Vk;
    Ilog(k) = Ik;
    Plog(k) = Vk * Ik;

    % 2) 在当前 D 基础上加一个很小的探测扰动，用于估算 dI/dV
    D_probe = min(max(D + test_step, Dmin), Dmax);
    [V_probe, I_probe] = solve_pv_operating_point(D_probe, Rload, G, Tc, pv);

    dV = V_probe - Vk;
    dI = I_probe - Ik;

    % 默认不调整占空比
    dD = 0;

    if abs(Vk) < epsV
        % 电压过小，避免除零，此时不更新 D
        dD = 0;
    else
        if abs(dV) < epsV
            % 近似 dV ≈ 0 的情况
            if abs(dI) < epsI
                % dV≈0 且 dI≈0 -> 认为已经非常接近 MPP
                dD = 0;
            elseif dI > 0
                % 电压没怎么变，电流增加 -> 在 MPP 左侧，需要 V↑
                % 对 Boost: V↑ => D↓
                dD = -DeltaD;
            else
                % 电压没怎么变，电流减小 -> 在 MPP 右侧，需要 V↓
                % 对 Boost: V↓ => D↑
                dD = +DeltaD;
            end
        else
            % 正常情况：使用增量电导公式
            dIdV = dI / dV;
            cond = dIdV + Ik / Vk;  % 近似 dI/dV + I/V

            if abs(cond) < epsCond
                % 认为已在 MPP 附近
                dD = 0;
            elseif cond > 0
                % 左侧: dP/dV > 0 => 需要提高 Vpv => 对 Boost: 减小 D
                dD = -DeltaD;
            else
                % 右侧: dP/dV < 0 => 需要降低 Vpv => 对 Boost: 增大 D
                dD = +DeltaD;
            end
        end
    end

    % 3) 更新占空比并限幅
    D = D + dD;
    if D > Dmax
        D = Dmax;
    elseif D < Dmin
        D = Dmin;
    end

    % 4) 用新的 D 解新的工作点 (给下一次迭代用)
    [Vk, Ik] = solve_pv_operating_point(D, Rload, G, Tc, pv);
end

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
    q = 1.602176634e-19; 
    k_b = 1.380649e-23;
    T_K = T + 273.15; 
    Vt  = params.Ns * k_b * T_K / q;
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
end