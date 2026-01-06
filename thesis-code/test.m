%% LFMCW 色散诊断参数敏感度验证仿真 (v2.0 - 鲁棒性增强版)
% 功能：
% 1. 验证群时延对电子密度(一阶)与碰撞频率(二阶)的敏感度差异
% 2. 【新增】加入幅度衰减监测，划定“信号可观测域”
% 3. 【改进】提升数值微分精度，优化观测频点选择

clc; clear; close all;

%% 1. 物理常数与仿真参数设置
c = 2.99792458e8;       % 光速 (m/s)
e = 1.60217663e-19;     % 电子电荷 (C)
me = 9.10938356e-31;    % 电子质量 (kg)
eps0 = 8.85418781e-12;  % 真空介电常数 (F/m)

d = 0.3;                % 等离子体有效厚度 (m)
f_start = 20e9;         % 扫频起始频率 30 GHz
f_stop = 40e9;          % 扫频终止频率 40 GHz
N_points = 10000;       % 【改进】点数提升至1万，消除微分截断误差
f = linspace(f_start, f_stop, N_points); 
omega = 2 * pi * f;     

% 系统动态范围阈值 (假设雷达接收灵敏度限制)
Loss_Threshold_dB = -80; 

%% 2. 基准状态设定
ne_0 = 5e18;            % 基准电子密度 (m^-3)
% 计算基准截止频率
wp_0 = sqrt(ne_0 * e^2 / (eps0 * me));
fp_0 = wp_0 / (2*pi);
fprintf('------------------------------------------------------\n');
fprintf('基准等离子体截止频率 fp = %.2f GHz\n', fp_0/1e9);
fprintf('------------------------------------------------------\n');

%% 3. 实验一：电子密度 ne 的敏感度测试 (时延主导性)
% 保持 nue 为典型值，变化 ne
nue_typ = 1e9; % 1 GHz
ne_list = [0.95 * ne_0, ne_0, 1.05 * ne_0]; 
legends_ne = {};

figure('Name', '参数敏感度与可观测性分析', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

% --- 子图1: 电子密度对群时延的影响 ---
subplot(2,2,1); hold on; grid on;
for i = 1:length(ne_list)
    ne_curr = ne_list(i);
    wp_curr = sqrt(ne_curr * e^2 / (eps0 * me));
    
    [tau_g, ~] = calculate_complex_drude(omega, wp_curr, nue_typ, d, c);
    
    plot(f/1e9, tau_g*1e9, 'LineWidth', 2);
    legends_ne{i} = sprintf('n_e = %.2e (%.0f%%)', ne_curr, (i-2)*5);
end
title('A. 电子密度对群时延的决定性作用');
xlabel('探测频率 (GHz)'); ylabel('群时延 (ns)');
legend(legends_ne, 'Location', 'NorthEast');
xlim([f_start/1e9 f_stop/1e9]);

%% 4. 实验二：碰撞频率 nue 的敏感度与衰减测试 (核心改进)
% 变化 nue，观察时延微扰 + 幅度衰减
nue_list = [0.1e9, 5e9, 10e9]; % 0.1GHz(弱), 5GHz(强), 20GHz(极端)
legends_nue = {};
colors = lines(length(nue_list));

% --- 子图2: 碰撞频率对群时延的微扰 ---
subplot(2,2,2); hold on; grid on;
% --- 子图3: 碰撞频率对幅度衰减的影响 ---
subplot(2,2,3); hold on; grid on;
plot([f_start/1e9, f_stop/1e9], [Loss_Threshold_dB, Loss_Threshold_dB], 'r--', 'LineWidth', 1.5); % 阈值线

for i = 1:length(nue_list)
    nue_curr = nue_list(i);
    [tau_g, loss_dB] = calculate_complex_drude(omega, wp_0, nue_curr, d, c);
    
    % 绘制时延 (Subplot 2)
    subplot(2,2,2);
    plot(f/1e9, tau_g*1e9, 'LineWidth', 2, 'Color', colors(i,:));
    
    % 绘制衰减 (Subplot 3)
    subplot(2,2,3);
    plot(f/1e9, loss_dB, 'LineWidth', 2, 'Color', colors(i,:));
    
    legends_nue{i} = sprintf('\\nu_e = %.1f GHz', nue_curr/1e9);
end

% 完善子图2
subplot(2,2,2);
title('B. 碰撞频率对群时延的二阶微扰');
xlabel('探测频率 (GHz)'); ylabel('群时延 (ns)');
legend(legends_nue, 'Location', 'NorthEast');
xlim([f_start/1e9 f_stop/1e9]);

% 完善子图3
subplot(2,2,3);
title('C. 碰撞频率对信号幅度的决定性影响');
xlabel('探测频率 (GHz)'); ylabel('透射衰减 S21 (dB)');
% 修改后（正确）：
legend([{'系统灵敏度阈值'}, legends_nue], 'Location', 'SouthEast');
ylim([-100, 0]); 
grid minor;

%% 5. 实验三：二阶小量验证 (修正版：解决标量微分NaN问题)
% 选取两个观测点：
% Point A: 靠近截止频率 (High Dispersion, 22 GHz)
% Point B: 典型诊断频点 (Ka-Band, 35 GHz)
f_obs_list = [22e9, 35e9]; 
line_styles = {'o-', 's--'};

subplot(2,2,4); hold on; grid on;

nue_range = logspace(8, 10.5, 30); % 0.1 GHz 到 30 GHz
legends_err = {};

fprintf('实验三：二阶小量验证数据\n');
for k = 1:length(f_obs_list)
    f_target = f_obs_list(k);
    
    % --- 关键修正 ---
    % 构造一个微小的频率向量 [f-delta, f, f+delta]
    % 这样 gradient 函数才能计算出中间点的导数
    delta_f = f_target * 0.001; % 0.1% 的步长
    f_local = [f_target - delta_f, f_target, f_target + delta_f];
    omega_local = 2 * pi * f_local;
    
    % 计算该频点的无损耗时延真值 (取中间点 index=2)
    [tau_vec_ref, ~] = calculate_complex_drude(omega_local, wp_0, 0, d, c);
    tau_ref = tau_vec_ref(2); 
    
    err_list = [];
    delta_list = [];
    
    for i = 1:length(nue_range)
        nue_val = nue_range(i);
        
        % 计算有损耗时延 (取中间点 index=2)
        [tau_vec_val, ~] = calculate_complex_drude(omega_local, wp_0, nue_val, d, c);
        tau_val = tau_vec_val(2);
        
        % 相对误差
        rel_err = abs(tau_val - tau_ref) / tau_ref;
        err_list(end+1) = rel_err;
        
        % 二阶小量 delta = (nue/omega)^2
        % 注意使用中心频率计算 delta
        delta_list(end+1) = (nue_val / (2*pi*f_target))^2;
    end
    
    loglog(delta_list, err_list, line_styles{k}, 'LineWidth', 1.5, 'MarkerSize', 6);
    legends_err{k} = sprintf('f_{obs} = %.0f GHz (f_p/f = %.2f)', f_target/1e9, fp_0/f_target);
    
    % 输出验证信息
    % 使用 interp1 插值获取特定点的误差数值用于打印
    fprintf('  频点 %.0f GHz: 当 nue=5GHz 时, 相对误差 = %.2e %%\n', ...
        f_target/1e9, interp1(nue_range, err_list, 5e9)*100);
end

% 绘制理论斜率参考线 (y = x)
% 为了让参考线位置好看，取数据的中间段
x_ref = [1e-4, 1e-1];
% 让参考线通过数据的某一个点，确保平行且靠近
ref_anchor_y = err_list(end); 
ref_anchor_x = delta_list(end);
y_ref = x_ref * (ref_anchor_y / ref_anchor_x); 

plot(x_ref, y_ref, 'k:', 'LineWidth', 2);
legends_err{end+1} = '理论斜率 1:1 (Linear)';

title('D. 误差与微扰项 \delta=(\nu_e/\omega)^2 的对应关系');
xlabel('微扰项 \delta'); ylabel('群时延相对误差');
legend(legends_err, 'Location', 'NorthWest');

%% --- 核心函数：全复数 Drude 模型计算 ---
function [tau_g, loss_dB] = calculate_complex_drude(omega, wp, nue, d, c)
    % 输入: 角频率, 特征角频率, 碰撞频率, 厚度, 光速
    % 输出: 群时延(s), 幅度衰减(dB)
    
    % 1. 复介电常数 (Complex Permittivity)
    % eps = 1 - wp^2 / (omega * (omega - j*nue))
    denom = omega.^2 + nue.^2;
    eps_real = 1 - (wp^2 ./ denom);
    eps_imag = - (nue ./ omega) .* (wp^2 ./ denom); % 虚部为负表示损耗(工程符号)
    
    eps_complex = eps_real + 1j * eps_imag;
    
    % 2. 复传播常数 gamma = alpha + j*beta
    % gamma = j * (omega/c) * sqrt(eps)
    % 注意 MATLAB sqrt 对复数直接开方得到 a+bi, 对应 gamma/(j*omega/c)
    % 也就是 sqrt(eps) = (beta - j*alpha) * (c/omega) ... 需小心符号
    
    % 更稳健的方法：直接求复折射率 N = n - j*k
    N = sqrt(eps_complex); 
    n_ref = real(N);      % 折射率实部
    kappa = -imag(N);     % 消光系数 (注意 MATLAB sqrt 结果虚部为负时，kappa为正)
    
    % 3. 相位常数 beta 与 衰减常数 alpha
    beta = (omega ./ c) .* n_ref;
    alpha = (omega ./ c) .* kappa;
    
    % 4. 计算群时延 (数值微分)
    % 使用中心差分 gradient，点数够多时精度极高
    d_omega = gradient(omega);
    d_beta = gradient(beta);
    tau_g = d .* (d_beta ./ d_omega);
    
    % 5. 计算幅度衰减 (S21 in dB)
    % 功率传输系数 T = exp(-2 * alpha * d)
    % dB = 10 * log10(T) = -20 * alpha * d * log10(e)
    loss_dB = 20 * log10(exp(-alpha .* d));
end