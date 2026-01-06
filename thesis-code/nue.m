%% LFMCW 等离子体诊断论文 - 3.2.3 参数敏感性分析仿真 (最终修正版)
% 修复内容：
% 1. 移除了 ylim 强制限制，解决“数据在楼下，窗口在楼上”导致的空白图问题。
% 2. 修复了 legend_str 未定义变量错误。
% 3. 优化了坐标轴显示 (axis tight)。

clc; clear; close all;

%% 1. 全局物理常量与仿真参数
c = 2.99792458e8;           % 光速 (m/s)
e = 1.60217663e-19;         % 电子电荷 (C)
me = 9.10938356e-31;        % 电子质量 (kg)
eps0 = 8.85418781e-12;      % 真空介电常数 (F/m)

% 仿真频段 (Ka波段 34-38 GHz)
f_start = 20e9;
f_end = 40e9;
N_points = 1000;            
f = linspace(f_start, f_end, N_points);
omega = 2 * pi * f;

% 介质几何参数 (仅计算介质本身，不含空气隙)
d = 0.15;                   % 等离子体厚度 (m)

%% 2. 基础参数设定
fp_base = 29e9;             
wp_base = 2 * pi * fp_base;
ne_base = (wp_base^2 * eps0 * me) / e^2; 
nu_base = 1.5e9;            

%% =============================================================
%% 3. 仿真 A: 电子密度 ne 的主导性分析 (Figure 3-3a)
%% =============================================================
figure('Name', 'Fig 3-3a: Electron Density Sensitivity', 'Color', 'w', 'Position', [100, 100, 600, 450]);

% 设定 ne 变化范围
ne_scales = [0.9, 1.0, 1.1];
colors = {'b', 'k', 'r'}; 
line_styles = {'--', '-', '-.'};

hold on; grid on;
% 【关键修正】: 必须先初始化图例容器！
legend_str = {}; 

for i = 1:length(ne_scales)
    ne_current = ne_base * ne_scales(i);
    
    % 计算
    [tau_g, ~] = calculate_drude_response(omega, ne_current, nu_base, d, c, eps0, me, e);
    
    % 绘图 (频率单位 GHz, 时延单位 ns)
    plot(f/1e9, tau_g*1e9, 'Color', colors{i}, 'LineWidth', 2, 'LineStyle', line_styles{i});
    
    % 记录图例
    fp_current = sqrt(ne_current * e^2 / (eps0 * me)) / (2*pi);
    legend_str{end+1} = sprintf('n_e (对应 f_p = %.1f GHz)', fp_current/1e9);
end

xlabel('探测频率 (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('群时延 (ns)', 'FontSize', 12, 'FontWeight', 'bold');
title('图 3-3a: 电子密度对群时延曲线的拓扑控制', 'FontSize', 13);

% 【关键修正】: 自动适应数据范围
axis tight; 
% 稍微留点边距美观一点
yl = ylim; ylim([yl(1)*0.95, yl(2)*1.05]); 

legend(legend_str, 'Location', 'Best', 'FontSize', 11);


%% =============================================================
%% 4. 仿真 B: 碰撞频率 nu 的解耦特性分析 (Figure 3-3b)
%% =============================================================
figure('Name', 'Fig 3-3b: Collision Frequency Sensitivity', 'Color', 'w', 'Position', [750, 100, 600, 450]);

% 设定 nu 变化范围
nu_list = [0.1e9, 1.5e9, 3e9];
colors_nu = {'[0 0.6 0]', 'k', 'm'}; 

% 【关键修正】: 初始化第二个图例容器！
legend_str_nu = {};

for i = 1:length(nu_list)
    nu_current = nu_list(i);
    
    [tau_g, mag_dB] = calculate_drude_response(omega, ne_base, nu_current, d, c, eps0, me, e);
    
    % 左轴画时延 (实线)
    yyaxis left;
    hold on;
    plot(f/1e9, tau_g*1e9, 'Color', colors_nu{i}, 'LineWidth', 2, 'LineStyle', '-');
    
    % 右轴画幅度 (虚线)
    yyaxis right;
    hold on;
    plot(f/1e9, mag_dB, 'Color', colors_nu{i}, 'LineWidth', 1.5, 'LineStyle', '--');
    
    legend_str_nu{end+1} = sprintf('\\nu_e = %.1f GHz', nu_current/1e9);
end

% 设置左轴属性
yyaxis left;
ylabel('群时延 (ns)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'k');
axis tight; % 让左轴自动适应数据 (约0.5-1.5ns)
yl_left = ylim; ylim([yl_left(1)*0.9, yl_left(2)*1.1]); % 微调

% 设置右轴属性
yyaxis right;
ylabel('透射幅度 S_{21} (dB)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'r');
ylim([-35, 5]); % 幅度范围通常比较固定

grid on;
xlabel('探测频率 (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
title('图 3-3b: 碰撞频率的“时延钝感”与“幅度敏感”解耦特性', 'FontSize', 13);

% 添加辅助说明
text(min(f)/1e9+0.2, min(yl_left)+0.1, '实线: 群时延 (左轴)', 'Parent', gca, 'FontSize', 10, 'Color', 'k');
% 注意：text默认坐标基于当前激活的轴(左轴)，如果想标右轴数据得小心，简单起见用legend区分颜色即可
legend(legend_str_nu, 'Location', 'SouthWest', 'FontSize', 11);


%% =============================================================
%% 5. 物理内核函数 (Drude 模型解析计算)
%% =============================================================
function [tau_g, mag_dB] = calculate_drude_response(omega, ne, nu, d, c, eps0, me, e)
    wp = sqrt(ne * e^2 / (eps0 * me));
    
    % 计算复介电常数
    term_denom = omega.^2 + nu^2;
    eps_real = 1 - (wp^2) ./ term_denom;
    eps_imag = -(nu ./ omega) .* (wp^2 ./ term_denom); 
    eps_complex = eps_real + 1i * eps_imag;
    
    % 计算复传播常数
    k0 = omega ./ c;
    gamma = 1i .* k0 .* sqrt(eps_complex);
    
    % 传递函数
    H = exp(-gamma * d);
    
    % 幅度 (dB)
    mag_dB = 20 * log10(abs(H));
    
    % 群时延 (数值微分)
    phi = unwrap(angle(H));
    d_phi = diff(phi);
    d_omega = diff(omega);
    tau_g_raw = -d_phi ./ d_omega;
    tau_g = [tau_g_raw, tau_g_raw(end)];
end