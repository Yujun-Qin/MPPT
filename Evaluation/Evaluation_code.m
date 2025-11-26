clear;clc;
file = dir('mppt_*.m');      % 找到算法文件
n = numel(file);             % 文件数目
mppt_func = cell(n,1);
name = cell(n,1);
E = zeros(n,1);              % 近似量化输出能量
P_part = cell(n,1);
tol = 3e-3;                  % 收敛容差

for i = 1:n
    [~,name{i},~] = fileparts(file(i).name);   % 提取文件名称
    mppt_func{i} = str2func(name{i});          % 将文件名称作为函数句柄
    
    t0 = cputime;                              % 计时CPU的运行时间
    f = mppt_func{i};
    P = f();                                   % 提取某算法对应的功率迭代过程
    dt = 1:length(P);
    % 判断迭代次数是否符合要求
    if length(dt) > 1000
        fprintf('迭代次数不符合要求')
        continue;
    end
    % 判断达到MPP的快慢以及是否收敛
    P_mpp = max(P);
    err = abs(P - P_mpp);                      
    tailMax = cummax(err(end:-1:1));             % 从后往前统计实时功率与最大功率的差值
    tailMax = tailMax(end:-1:1);                 % 再反转回来，得到从前往后的差值
    P_start = find(tailMax <= tol, 1, 'first');  % 找到最开始到达最大功率的点

    if isempty(P_start)
        fprintf('%s 未收敛到规定范围\n', name{i});
        continue;                                % 此算法不符合要求，处理下一个
    end

    ratio = P_start / length(P);                 % 计算收敛比例，越小越快   
    E(i) = trapz(dt, P);                     % 计算总能量
    t = timeit(@()f());
    T = cputime - t0;
    fprintf(['\n--- %s ---\n', ...
            '算法输出能量：%.6g\n',...
            '平均耗时: %.6g 秒\n',...
            '代码CPU耗时：%.4f 秒\n',...
            '收敛比例：%.4g \n\n'],...
            name{i}, E(i), t, T,ratio);
end