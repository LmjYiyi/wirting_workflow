% =========================================================================
%               LFMCW 色散效应工程判据可视化分析 (终极修正版)
%               对应论文章节：3.4.2 & 3.1.3
% =========================================================================

%% 0. 初始化
clc; clear; close all;

%% 1. 定义基本物理和仿真参数
c = 3e8;                  % 光速 (m/s)
d = 150e-3;               % 等离子体厚度 (m)
tau0 = d / c;             % 真空基准时延 (s)

% --- 定义用于总体分析的参数范围 ---
fp_values_general = [10e9, 20e9, 30e9]; % 截止频率
B_values_general  = [0.5e9, 1.0e9, 2.0e9]; % 雷达带宽

%% 2. 总体分析：色散判据 \xi vs 测量中心频率
fprintf('--- 正在生成总体分析图 ---\n');

figure('Name', '总体分析', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; grid on;
colors = lines(length(fp_values_general));
line_styles = {'-', '--', ':'};

for b_idx = 1:length(B_values_general)
    B = B_values_general(b_idx);
    
    for fp_idx = 1:length(fp_values_general)
        fp = fp_values_general(fp_idx);
        
        % 扫描频率范围：从 1.05倍截止频率 到 50 GHz
        f0_vec = linspace(fp * 1.02, 50e9, 1000); 
        
        % === 核心计算 (调用优化后的函数) ===
        Xi_vec = calculate_criterion_xi(f0_vec, fp, B, tau0);
        
        % 绘图
        plot(f0_vec/1e9, Xi_vec, ...
            'Color', colors(fp_idx, :), ...
            'LineStyle', line_styles{b_idx}, ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('f_p=%dG, B=%.1fG', round(fp/1e9), B/1e9));
    end
end

% 阈值线
yline(1, 'r-', 'LineWidth', 2, 'DisplayName', '工程阈值 \xi = 1');

% 区域标注
text(45, 0.5, '安全区 (Safe)', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
text(45, 2, '失效区 (Fail)', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

set(gca, 'YScale', 'log');
xlabel('测量中心频率 f_0 (GHz)');
ylabel('色散判据值 \xi = B \cdot \eta \cdot \tau_0');
title('总体分析: 不同截止频率与带宽下的色散强度分布');
legend('show', 'Location', 'northeast', 'NumColumns', 2);
ylim([1e-2, 1e2]);
xlim([10, 50]);

%% 3. 特定参数分析：Ka波段典型工况
fprintf('\n--- 正在生成特定参数分析图 ---\n');

% --- 典型工况参数 ---
fp_spec = 33.59e9;        % 高密度等离子体截止频率 (~1.4e19 m^-3)
B_spec_list = [0.5e9, 2.0e9, 4.0e9]; % 对比：窄带 vs 宽带 vs 超宽带
f0_start = 34e9;          % Ka波段起始
f0_stop  = 40e9;          % Ka波段终止

figure('Name', '特定参数分析', 'Color', 'w', 'Position', [950, 100, 600, 500]);
hold on; grid on;
colors_spec = {'b', 'k', 'r'};

f0_vec = linspace(f0_start, f0_stop, 500);

for i = 1:length(B_spec_list)
    B = B_spec_list(i);
    Xi_vec = calculate_criterion_xi(f0_vec, fp_spec, B, tau0);
    
    plot(f0_vec/1e9, Xi_vec, ...
        'Color', colors_spec{i}, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('带宽 B = %.1f GHz', B/1e9));
end

yline(1, 'm--', 'LineWidth', 2, 'DisplayName', '阈值 \xi=1');

set(gca, 'YScale', 'log');
xlabel('测量中心频率 f_0 (GHz)');
ylabel('色散判据值 \xi');
title(['特定工况分析 (f_p = ' num2str(fp_spec/1e9) ' GHz)']);
legend('Location', 'Best');
ylim([1e-1, 1e2]);

% 添加背景色块表示 Ka 波段常用区
x_patch = [34, 36, 36, 34];
y_patch = [1e-1, 1e-1, 1e2, 1e2];
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Ka-Band (34-36G)');

%% 4. 本地函数：严格对应论文公式 (3-18) 和 (3.4.2)
function Xi = calculate_criterion_xi(f, fp, B, tau0)
    % 输入: 
    % f  - 探测频率向量 (Hz)
    % fp - 截止频率 (Hz)
    % B  - 信号带宽 (Hz)
    % tau0 - 真空基准时延 (s)
    
    % 1. 计算归一化频率比 x = fp / f
    x = fp ./ f;
    
    % 2. 避免复数 (只计算 f > fp 的部分)
    valid_idx = f > fp;
    Xi = NaN(size(f)); % 截止区内设为 NaN
    
    if any(valid_idx)
        x_valid = x(valid_idx);
        f_valid = f(valid_idx);
        
        % 3. 计算非线性度因子 eta (对应论文式 3-18)
        % eta = (B/f) * [x^2 / (1-x^2)^1.5]
        term1 = B ./ f_valid;
        term2 = (x_valid.^2) ./ ((1 - x_valid.^2).^1.5);
        eta = term1 .* term2;
        
        % 4. 计算工程判据 Xi (对应论文 3.4.2)
        % Xi = B * eta * tau0
        Xi(valid_idx) = B .* eta .* tau0;
    end
end