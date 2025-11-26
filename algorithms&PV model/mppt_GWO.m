function [Plog] = mppt_GWO()

% ==== MPPT：灰狼优化算法（GWO-MPPT）====

Nstep   = 100;    % 迭代代数
Nw      = 8;      % 狼的数量
Dlog    = zeros(1,Nstep);
Vlog    = zeros(1,Nstep);
Ilog    = zeros(1,Nstep);
Plog    = zeros(1,Nstep);

% 初始化狼群位置
pos   = Dmin + (Dmax - Dmin) * rand(Nw,1); 
fit   = -1e9 * ones(Nw,1);

% 初始化 α, β, δ
alpha_pos = pos(1); alpha_fit = -1e9;
beta_pos  = pos(1); beta_fit  = -1e9;
delta_pos = pos(1); delta_fit = -1e9;

for k = 1:Nstep
    % 1) 评估所有狼当前适应度
    for i = 1:Nw
        Di = pos(i);
        [Vi, Ii] = solve_pv_operating_point(Di, Rload, G, Tc, pv);
        Pi = Vi * Ii;
        fit(i) = Pi;

        % 更新 α, β, δ
        if Pi > alpha_fit
            delta_fit = beta_fit; delta_pos = beta_pos;
            beta_fit  = alpha_fit; beta_pos  = alpha_pos;
            alpha_fit = Pi;        alpha_pos = Di;
        elseif Pi > beta_fit
            delta_fit = beta_fit; delta_pos = beta_pos;
            beta_fit  = Pi;       beta_pos  = Di;
        elseif Pi > delta_fit
            delta_fit = Pi;       delta_pos = Di;
        end
    end

    % 2) 用 α 对应的占空比作为本代输出并记录
    D  = alpha_pos;
    [Vk, Ik] = solve_pv_operating_point(D, Rload, G, Tc, pv);
    Pk       = Vk * Ik;
    Dlog(k)  = D;
    Vlog(k)  = Vk;
    Ilog(k)  = Ik;
    Plog(k)  = Pk;

    % 3) 更新参数 a (线性递减从 2 -> 0)
    a = 2 * (1 - k / Nstep);
    if a < 0, a = 0; end

    % 4) 根据 α, β, δ 更新所有狼的位置
    for i = 1:Nw
        X = pos(i);

        r1 = rand(); r2 = rand();
        A1 = 2*a*r1 - a;
        C1 = 2*r2;
        D_alpha = abs(C1 * alpha_pos - X);
        X1 = alpha_pos - A1 * D_alpha;

        r1 = rand(); r2 = rand();
        A2 = 2*a*r1 - a;
        C2 = 2*r2;
        D_beta = abs(C2 * beta_pos - X);
        X2 = beta_pos - A2 * D_beta;

        r1 = rand(); r2 = rand();
        A3 = 2*a*r1 - a;
        C3 = 2*r2;
        D_delta = abs(C3 * delta_pos - X);
        X3 = delta_pos - A3 * D_delta;

        X_new = (X1 + X2 + X3) / 3;

        % 限幅到 [Dmin, Dmax]
        if X_new > Dmax
            X_new = Dmax;
        elseif X_new < Dmin
            X_new = Dmin;
        end

        pos(i) = X_new;
    end
end


% ==== 迭代过程曲线 ====
figure; 
plot(Plog,'LineWidth',1.2); grid on; xlabel('Iteration'); ylabel('P (W)');
title('Iteration of Power(GWO)');

figure; 
plot(Vlog,'LineWidth',1.2); grid on; xlabel('Iteration'); ylabel('V_{pv} (V)');
title('Iteration of Voltoge(GWO)');

figure; 
plot(Dlog,'LineWidth',1.2); grid on; xlabel('Iteration'); ylabel('D');
title('Iteration of Duty(GWO)');

end