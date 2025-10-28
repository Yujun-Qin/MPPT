function [V, I, P] = pv_model(Vpv, scenario, time, pv_params)
% 输入: Vpv - 光伏电压, scenario - 测试场景, time - 当前时间
% 输出: V, I, P - 电压,电流,功率
persistent G T;   
if isempty(G)
   G = 1000; % 标准光照(W/m²)
   T = 25;   % 标准温度(°C)
end

% 根据场景调整环境参数
switch scenario
   case 1 % 恒定光照
   G_actual = G;
   T_actual = T;           
   case 2 % 光照突变
        if time < 1
           G_actual = 1000;
        else
           G_actual = 600; % 光照突变
        end
   T_actual = T;          
   case 3 % 局部阴影 - 双峰特性，模拟局部阴影下的多峰P-V曲线
   Vmp = pv_params.Vmp;
   if Vpv < Vmp * 0.7
      G_actual = 800;
   else
      G_actual = 400; % 阴影区域
   end
      T_actual = T;         
   case 4 % 温度变化
       G_actual = G;
       T_actual = 25 + 20 * sin(2*pi*0.5*time); % 温度波动         
   otherwise
       G_actual = G;
       T_actual = T;
end
    
    % 简化光伏模型 (可根据需要替换为更精确的模型)
    Isc = pv_params.Isc * (G_actual/1000);
    Voc = pv_params.Voc + 0.05*(T_actual-25);
    
    % 单二极管模型近似
    a = 1.5; % 理想因子
    Vt = 1.38e-23 * (273+T_actual) / 1.6e-19 * pv_params.Ns;
    
    I = Isc * (1 - exp((Vpv - Voc)/(a*Vt)));
    
    % 对于局部阴影场景，创建双峰特性
    if scenario == 3
        I2 = 0.6 * Isc * (1 - exp((Vpv - 0.6*Voc)/(a*Vt)));% 添加第二个峰值
        I = max(I, I2);
    end   
    P = Vpv .* I;
    
    % 添加测量噪声
    I = I + 0.01*Isc*randn(size(I));
end