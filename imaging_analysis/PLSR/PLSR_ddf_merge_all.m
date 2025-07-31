clear; clc; close all;
%% Base parameters
y_limit = [-1.5, 1.5];
data_range = [-3, 3];
fix_node_size = 00; % fix atlas plot node size. if not to use, set this value 0. recommend 500
sym_option = true;
alpha = 0.1;
zscore_option = true;
intervals = {'start', 'continue', 'stop'}; interval_all = [{'resting'}, intervals];

figure_plot_option = true; % plot node-wise figure;
%% Choose task
addpath(genpath('E:\Final Codes'))
task_list = {'treadmill', 'wheel', 'disk'};
task_choice = questdlg('Select Task:', 'Task Selection', ...
    'treadmill', 'wheel', 'disk', 'treadmill'); % default == 'treadmill'
if isempty(task_choice)
    disp('Task selection canceled. Exiting...');
    return;
end
task_map = containers.Map({'treadmill', 'wheel', 'disk'}, {'treadmill', 'wheel', 'disk'});
task = task_map(task_choice);

disp(['Selected Task: ', task]);

resultDir = fullfile('D:\data analysis\results', 'PLSR_ddf_all');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
taskDir = fullfile(resultDir, task);
if ~exist(taskDir, 'dir')
    mkdir(taskDir);
end
colormap_list = load('D:\reference\colormap_data.mat');
%% Load files
parentFolder = 'D:\data analysis\locomotion training';
folderList = dir(parentFolder);
folderList = folderList([folderList.isdir]);
folderList = folderList(~ismember({folderList.name}, {'.', '..'}));

[idx, tf] = listdlg('PromptString', 'Select folders:',...
                   'ListString', {folderList.name},...
                   'SelectionMode', 'multiple');

if ~tf
    return;
end

topFolders = fullfile(parentFolder, {folderList(idx).name});

numFolders = numel(topFolders);
data = struct('mouseName', {}, 'data', {});
avg_dzs_1x = cell(numFolders, 1);
avg_dzs_2x = cell(numFolders, 1);
resting_1x = cell(numFolders, 1);
resting_2x = cell(numFolders, 1);
avg_zs_1x = cell(numFolders, 1);
avg_zs_2x = cell(numFolders, 1);
zscore_1x = cell(numFolders, 1);
zscore_2x = cell(numFolders, 1);

for i = 1:numFolders
    topFolder = topFolders{i};
    [~, mouseName] = fileparts(topFolder);
    subFolders = dir(topFolder);
    subFolders = subFolders([subFolders.isdir]);
    subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

    [~, idx_recent] = max([subFolders.datenum]);
    recentFolder = fullfile(topFolder, subFolders(idx_recent).name);
    if zscore_option == true
        targetFile = fullfile(recentFolder, 'imaging', 'df_data','PLSR_all','PLSR_ddf_zscored.mat');
        subfix = '_zscored';
    elseif zscore_option == false
        targetFile = fullfile(recentFolder, 'imaging', 'df_data','PLSR_all','PLSR_ddf.mat');
        subfix = '';
    end
    if exist(targetFile, 'file')
        loadedData = load(targetFile);
        data(i).mouseName = mouseName;
        data(i).data = loadedData;

        if isfield(loadedData, 'mean_dzs_1x')
            avg_dzs_1x{i} = loadedData.mean_dzs_1x; 
            avg_dzs_2x{i} = loadedData.mean_dzs_2x; 
            resting_1x{i} = loadedData.zs_rest_1x; 
            resting_2x{i} = loadedData.zs_rest_2x;
            avg_zs_1x{i} = loadedData.mean_zs_1x; 
            avg_zs_2x{i} = loadedData.mean_zs_2x; 
            zscore_1x{i} = loadedData.zscore_1x; 
            zscore_2x{i} = loadedData.zscore_2x; 
        else
            warning('data not found in file %s', targetFile);
        end
    else
        warning('File does not exist: %s', targetFile);
    end
end
save(fullfile(taskDir, 'raw_data.mat'), "avg_dzs_1x", "avg_dzs_2x", "avg_zs_1x", "avg_zs_2x", "resting_1x", "resting_2x");
%% Calculate mean and std for each dataset
ddf_mean = mean_cell_matrices(avg_dzs_1x); ddf_mean2 = mean_cell_matrices(avg_dzs_2x);
df_mean = mean_cell_matrices(avg_zs_1x); df_mean2 = mean_cell_matrices(avg_zs_2x);
resting_mean = mean_cell_matrices(resting_1x); resting_mean2 = mean_cell_matrices(resting_2x);
ddf_std = std_cell_matrices(avg_dzs_1x); ddf_se = sem_cell_matrices(avg_dzs_1x); 

load('D:\reference\sorted_node_positions.mat');
load('D:\reference\reference_atlas.mat');
save(fullfile(taskDir, 'mean_data.mat'), 'ddf_mean','ddf_mean2','df_mean','df_mean2','resting_mean','resting_mean2', 'ddf_std', 'ddf_se');
%% ΔF/F 값에 대한 t-test 수행 (resting 기준)
% ΔF = df_mean - resting_mean (animal 단위)

% 가정: avg_zs_1x = {N x 1 cell} → 각 cell: [interval x nodes]
%       resting_1x = {N x 1 cell} → 각 cell: [1 x nodes]
% 기본 설정
num_animals = numel(avg_zs_1x);
n_ints = size(avg_zs_1x{1}, 1);
n_units = size(avg_zs_1x{1}, 2);

%% Filt data t-test
% 1x trial t-test
[p_valuesF1, h_valuesF1, sorted_pF1, sorted_idxF1] = run_ttest_fdr(avg_dzs_1x, alpha);

% 2x trial t-test
[p_valuesF2, h_valuesF2, sorted_pF2, sorted_idxF2] = run_ttest_fdr(avg_dzs_2x, alpha);
%% save data
save(fullfile(taskDir, sprintf('ttest_%s.mat', task)), ...
    'p_valuesF1', 'h_valuesF1', 'sorted_idxF1', 'sorted_pF1', ...
    'p_valuesF2', 'h_valuesF2', 'sorted_idxF2', 'sorted_pF2');
% save(fullfile(taskDir, sprintf(['ttest_%s', subfix,'.mat'], task)), ...
%     'p_values_zs1', 'h_valuesF1', 'sorted_idx_zs1', 'sorted_p_values_zs1', 'crit_p_zs1');
% return;

%% Plot averaged data: same range
mean_org1x = mean_cell_matrices(zscore_1x);
% 1x trial img for unit-wise
figure;
set(gcf, 'WindowState','maximized');
org_data_1x = mean_org1x(11:370, :);

imagesc(org_data_1x')
title('Original Data (1x)')
colorbar;
colormap('jet');
caxis([data_range(1), data_range(2)])
xline(80.5, '--k', 'LineWidth', 1.5);
xline(280.5, '--k', 'LineWidth', 1.5);
draw_region_boundaries(size(org_data_1x, 2));

saveas(gcf, fullfile(taskDir, 'PLSR_results_1x_heatmap.png'));
saveas(gcf, fullfile(taskDir, 'PLSR_results_1x_heatmap.fig'));
plot_epsname = fullfile(taskDir, 'PLSR_results_1x_heatmap.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
% return;
%% 
raw_animals1 = cell(num_animals, 1); masked_animals1 = cell(num_animals, 1);
for aidx = 1:num_animals
    curr1_mean_raw = NaN(n_ints, n_units/2); 
    curr1_mean_mask = NaN(n_ints, n_units/2);
    curr1 = avg_dzs_1x{aidx}; curr1_org = curr1; curr1(isnan(h_valuesF1')) = NaN;

    for eidx = 1:n_ints
        curr1_mean_raw(eidx, :) = average_pairs(curr1_org(eidx,:));
        curr1_mean_mask(eidx, :) = average_pairs(curr1(eidx,:));
    end
    raw_animals1{aidx} = curr1_mean_raw; masked_animals1{aidx} = curr1_mean_mask;
end
save(fullfile(taskDir, sprintf('PLSR_alldf_%s%s.mat', task, subfix)), 'raw_animals1', 'masked_animals1');
% winopen(taskDir)
%% scatter plot
data = ddf_mean;
ordered_regions = {'M2', 'M1', 'S1', 'Aud', 'Vis', 'RSC'};
num_regions = numel(ordered_regions);
colors = lines(num_regions);
int_name = {'Start', 'Continue', 'Stop'};  % 시간 포인트 이름
axis_limits = [
    0.0, 1.5,  0.0, 1.5;   % t = 1 (Start)
    0.0, 1.0,  0.0, 1.0;   % t = 2 (Continue)
    -1.0, 0.0,  -1.0, 0.0    % t = 3 (Stop)
];

% region_labels 생성
region_labels = strings(1, 28);
region_breaks = [6, 10, 16, 18, 24];
region_names = ordered_regions;
start_idx = 1;
for i = 1:length(region_breaks)
    stop_idx = region_breaks(i);
    region_labels(start_idx:stop_idx) = region_names{i};
    start_idx = stop_idx + 1;
end
region_labels(start_idx:end) = region_names{end};  % 마지막 RSC

% 반복문: 시간 포인트별로 그리기 및 저장
for t = 1:3
    % 좌우 데이터 추출
    left_idx = 1:2:27;
    right_idx = 2:2:28;
    x = data(t, left_idx);
    y = data(t, right_idx);
    regions = region_labels(left_idx);

    % 영역 인덱싱
    region_idx = zeros(size(regions));
    for i = 1:numel(ordered_regions)
        region_idx(strcmp(regions, ordered_regions{i})) = i;
    end

    % Figure 생성
    fig = figure('Visible', 'on'); hold on;
    for i = 1:length(x)
        c = colors(region_idx(i), :);
        scatter(x(i), y(i), 80, 'filled', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', c);
    end

    % 범례 생성
    h = gobjects(num_regions, 1);
    for i = 1:num_regions
        h(i) = scatter(nan, nan, 80, 'filled', ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'k');
    end
    legend(h, ordered_regions, 'Location', 'bestoutside');

    % 라벨 및 축
    xlabel('Left hemisphere');
    ylabel('Right hemisphere');
    title(sprintf('%s period', int_name{t}));
    axis equal;
    
    % 자동 축 범위 (0.5 단위 맞춤)
    min_val = min([x, y]);
    max_val = max([x, y]);
    min_lim = floor(min_val / 0.5) * 0.5;
    max_lim = ceil(max_val / 0.5) * 0.5;
    xlim([min_lim, max_lim]);
    ylim([min_lim, max_lim]);
    
    % y = x 기준선
    xline = linspace(min_lim, max_lim, 100);
    plot(xline, xline, 'k--', 'LineWidth', 1.2);

    % 회귀 분석
    mdl = fitlm(x', y');  % row vector -> column vector
    
    % 회귀선 시각화
    x_fit = linspace(min_lim, max_lim, 100);
    y_fit = predict(mdl, x_fit');
    plot(x_fit, y_fit, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);  % 회색 실선
    
    % 상관계수 계산
    R = corrcoef(x, y);  % 2x2 행렬
    r_value = R(1, 2);
    
    % R² 표시
    text_x = min_lim + 0.05;
    text_y = max_lim - 0.05;
    text(text_x, text_y, sprintf('R² = %.3f', mdl.Rsquared.Ordinary), ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    % (옵션) r-value도 함께 표시
    text(text_x, text_y - 0.05, sprintf('r = %.3f', r_value), ...
        'FontSize', 10, 'FontWeight', 'normal');

    grid on;

    % 저장
    fname = sprintf('Scatter_LeftRight_%s', int_name{t});
    set(fig,'renderer','Painters');
    print('-depsc','-tiff','-r300','-painters', fullfile(taskDir, [fname '.eps']));
    print(fig, fullfile(taskDir, [fname '.png']), '-dpng', '-r300');
    saveas(fig, fullfile(taskDir, [fname '.fig']));
    close(fig);
end
%% 
for t = 1:3
    fig = figure('Visible', 'on'); hold on;

    % 전체 좌/우 데이터 누적용
    all_x = [];
    all_y = [];
    all_region_idx = [];

    % 반복: 각 마우스
    for m = 1:5
        data_mouse = avg_dzs_1x{m};  % [3 x 28]
        x = data_mouse(t, 1:2:27);  % 좌
        y = data_mouse(t, 2:2:28);  % 우

        % region 인덱싱
        regions = region_labels(1:2:27);
        region_idx = zeros(size(regions));
        for i = 1:numel(ordered_regions)
            region_idx(strcmp(regions, ordered_regions{i})) = i;
        end

        % 점 추가
        for i = 1:length(x)
            c = colors(region_idx(i), :);
            scatter(x(i), y(i), 100, 'filled', ...
                'MarkerFaceColor', c, ...
                'MarkerFaceAlpha', 1.0, ...
                'MarkerEdgeColor', 'none');  % ✅ 외곽선 제거
        end

        all_x = [all_x, x];
        all_y = [all_y, y];
        all_region_idx = [all_region_idx, region_idx];
    end

    % 범례 생성
    h = gobjects(num_regions, 1);
    for i = 1:num_regions
        h(i) = scatter(nan, nan, 100, 'filled', ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'none');  % ✅ 범례 점도 외곽선 제거
    end
    legend(h, ordered_regions, 'Location', 'bestoutside');

    % 축 설정 (0.5 단위로 정리)
%     min_val = min([all_x, all_y]);
%     max_val = max([all_x, all_y]);
%     min_lim = floor(min_val / 0.5) * 0.5;
%     max_lim = ceil(max_val / 0.5) * 0.5;
    min_val = [-0.5, -0.5, -2.0];
    max_val = [2.0, 2.1, 0.5];
    min_lim = min_val(t); max_lim = max_val(t);

    xlim([min_lim, max_lim]);
    ylim([min_lim, max_lim]);
    axis equal;
    grid on;

    % y = x 기준선
    xx = linspace(min_lim, max_lim, 100);
    plot(xx, xx, 'k--', 'LineWidth', 1.2);

    % 회귀 분석
    mdl = fitlm(all_x', all_y');
    x_fit = linspace(min_lim, max_lim, 100);
    y_fit = predict(mdl, x_fit');
    plot(x_fit, y_fit, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

    % 상관계수 계산
    R = corrcoef(all_x, all_y);
    r_value = R(1, 2);

    % 텍스트 위치
    text_x = min_lim + 0.05;
    text_y = max_lim - 0.05;

    % R² 및 r 표시
    text(text_x, text_y, sprintf('R² = %.3f', mdl.Rsquared.Ordinary), ...
        'FontSize', 10, 'FontWeight', 'bold');
    text(text_x, text_y - 0.1, sprintf('r = %.3f', r_value), ...
        'FontSize', 10, 'FontWeight', 'normal');

    xlabel('Left hemisphere');
    ylabel('Right hemisphere');
    title(sprintf('All animals (Δz-score, %s)', int_name{t}));

    % 저장
    fname = sprintf('Scatter_AllAnimals_UnitPair_%s', int_name{t});
    set(fig,'renderer','Painters');
    print('-depsc','-tiff','-r300','-painters', fullfile(taskDir, [fname '.eps']));
    print(fig, fullfile(taskDir, [fname '.png']), '-dpng', '-r300');
    saveas(fig, fullfile(taskDir, [fname '.fig']));
    close(fig);
end
%% 
save_dir = taskDir;  % 저장할 디렉토리

for t = 1:3
    interval_name = int_name{t};  % 'start', 'continue', 'stop'
    
    % track 이름 (필요시 확장 가능)
    track_name = task;  % 예: 1x trial 기준

    % 조건별 데이터 초기화
    Subject = []; Left = []; Right = [];

    for m = 1:5
        data_mouse = avg_dzs_1x{m};  % [3 x 28]
        x = data_mouse(t, 1:2:27);  % 좌
        y = data_mouse(t, 2:2:28);  % 우

        n = numel(x);
        Subject = [Subject; (1:n)'];
        Left = [Left; x(:)];
        Right = [Right; y(:)];
    end

    % 테이블 생성
    paired_table = table(Subject, Left, Right);

    % 저장: .mat
    fname = sprintf('prism_%s_%s.mat', track_name, interval_name);
    save(fullfile(save_dir, fname), 'paired_table');

    % 저장: .xlsx (Prism에 바로 붙여넣기 용이)
    writetable(paired_table, fullfile(save_dir, [fname(1:end-4), '.xlsx']));
end

return;
%% Draw Bar plot
left_idx = 1:2:27;
right_idx = 2:2:28;
num_pairs = length(left_idx);
int_name = {'Start', 'Continue', 'Stop'};
alpha = 0.05;  % 유의수준

for t = 1:3
    % 평균과 표준편차
    x_mean = ddf_mean(t, left_idx);
    y_mean = ddf_mean(t, right_idx);
    x_std  = ddf_se(t, left_idx);
    y_std  = ddf_se(t, right_idx);

    % 평균과 std 배열 정리
    means = [x_mean; y_mean];  % [2 x num_pairs]
    stds  = [x_std;  y_std];

    pvals_raw = zeros(1, num_pairs);
    for i = 1:num_pairs
        left_unit  = left_idx(i);
        right_unit = right_idx(i);
    
        left_vals  = zeros(5,1);
        right_vals = zeros(5,1);
        for m = 1:5
            data_mouse = avg_dzs_1x{m};  % [3 x 28]
            left_vals(m)  = data_mouse(t, left_unit);
            right_vals(m) = data_mouse(t, right_unit);
        end
    
        [~, p] = ttest(left_vals, right_vals);  % paired t-test
        pvals_raw(i) = p;
    end
    
    % ✅ FDR 보정 (Benjamini-Hochberg)
    [p_sorted, sort_idx] = sort(pvals_raw);
    n = length(pvals_raw);
    fdr_thresh = (1:n) * alpha / n;
    below_thresh = find(p_sorted <= fdr_thresh, 1, 'last');
    
    % 보정된 유의성 표시
    sig_idx = false(1, n);
    if ~isempty(below_thresh)
        sig_idx(sort_idx(1:below_thresh)) = true;
    end

    % Figure
    fig = figure('Visible', 'on'); hold on;
    b = bar(means');  % [num_pairs x 2]
    b(1).FaceColor = [0.4 0.6 1];  % Left
    b(2).FaceColor = [1 0.5 0.5];  % Right

    % Error bars
    groupwidth = min(0.8, 2/(2 + 1.5));
    for i = 1:2
        x_pos = (1:num_pairs) - groupwidth/2 + (2*i - 1) * groupwidth / (2*2);
        errorbar(x_pos, means(i,:), stds(i,:), 'k', 'linestyle', 'none');
    end

    % 유의성 표시 (*, **, ***)
    for i = 1:num_pairs
        if sig_idx(i)
            p = pvals_raw(i);
            y_max = max(means(:,i) + stds(:,i)) + 0.05;
            x1 = i - groupwidth/4;
            x2 = i + groupwidth/4;
            plot([x1, x2], [y_max, y_max], 'k-', 'LineWidth', 1.2);
    
            % 별 개수 설정
            if p < 0.001
                stars = '***';
            elseif p < 0.01
                stars = '**';
            else
                stars = '*';
            end
            text(i, y_max + 0.02, stars, 'HorizontalAlignment', 'center', 'FontSize', 14);
        end
    end


    % X축 라벨
    xticks(1:num_pairs);
    xticklabels(arrayfun(@(l,r) sprintf('%d/%d', l, r), left_idx, right_idx, 'UniformOutput', false));
    xtickangle(0);

    legend({'Left', 'Right'});
    ylabel('Δz-score');
    title(sprintf('Left vs Right per unit (Δz-score, %s)', int_name{t}));

    y_max_overall = max(max(means + stds));
    y_min_overall = min(min(means - stds));
    y_min = floor(y_min_overall / 0.5) * 0.5;
    y_max = ceil(y_max_overall / 0.5) * 0.5;
    ylim([y_min, y_max]);
    yticks(y_min:0.5:y_max);

    grid on;

    % 저장
    fname = sprintf('Bar_Unitwise_LR_ddfmean_%s', int_name{t});
    set(fig,'renderer','Painters');
    print('-depsc','-tiff','-r300','-painters', fullfile(taskDir, [fname '.eps']));
    print(fig, fullfile(taskDir, [fname '.png']), '-dpng', '-r300');
    saveas(fig, fullfile(taskDir, [fname '.fig']));
    close(fig);
end
%% 
% 영역 정의
region_map = {
    'M2', 1:6;
    'M1', 7:10;
    'S1', 11:16;
    'Aud', 17:18;
    'Vis', 19:24;
    'RSC', 25:28;
};

ordered_regions = region_map(:,1);
num_regions = numel(ordered_regions);
alpha = 0.05;

for t = 1:3
    % 좌/우 영역 평균 계산
    left_means = zeros(1, num_regions);
    right_means = zeros(1, num_regions);
    left_se = zeros(1, num_regions);
    right_se = zeros(1, num_regions);
    pvals_raw = zeros(1, num_regions);

    for r = 1:num_regions
        units = region_map{r, 2};
        left_units = units(mod(units, 2) == 1);   % 홀수 → 좌
        right_units = units(mod(units, 2) == 0);  % 짝수 → 우

        left_vals = zeros(num_animals, 1);
        right_vals = zeros(num_animals, 1);

        for m = 1:num_animals
            data_mouse = avg_dzs_1x{m};  % [3 x 28]
            left_vals(m) = mean(data_mouse(t, left_units), 'omitnan');
            right_vals(m) = mean(data_mouse(t, right_units), 'omitnan');
        end

        left_means(r) = mean(left_vals, 'omitnan');
        right_means(r) = mean(right_vals, 'omitnan');
        left_se(r) = std(left_vals, 'omitnan') / sqrt(num_animals);
        right_se(r) = std(right_vals, 'omitnan') / sqrt(num_animals);

        [~, pvals_raw(r)] = ttest(left_vals, right_vals);
    end

    % ✅ FDR 보정
    [p_sorted, sort_idx] = sort(pvals_raw);
    fdr_thresh = (1:num_regions) * alpha / num_regions;
    below_thresh = find(p_sorted <= fdr_thresh, 1, 'last');

    sig_idx = false(1, num_regions);
    if ~isempty(below_thresh)
        sig_idx(sort_idx(1:below_thresh)) = true;
    end

    % Plot
    fig = figure('Visible', 'on'); hold on;
    means = [left_means; right_means];
    errs = [left_se; right_se];

    b = bar(means');  % [num_regions x 2]
    b(1).FaceColor = [0.4 0.6 1];  % Left
    b(2).FaceColor = [1 0.5 0.5];  % Right

    groupwidth = min(0.8, 2/(2 + 1.5));
    for i = 1:2
        x_pos = (1:num_regions) - groupwidth/2 + (2*i - 1) * groupwidth / (2*2);
        errorbar(x_pos, means(i,:), errs(i,:), 'k', 'linestyle', 'none');
    end

    % 유의성 (*, **, ***) 표시
    for i = 1:num_regions
        if sig_idx(i)
            p = pvals_raw(i);
            y_max = max(means(:,i) + errs(:,i)) + 0.05;
            x1 = i - groupwidth/4;
            x2 = i + groupwidth/4;
            plot([x1, x2], [y_max, y_max], 'k-', 'LineWidth', 1.2);

            if p < 0.001
                stars = '***';
            elseif p < 0.01
                stars = '**';
            else
                stars = '*';
            end
            text(i, y_max + 0.02, stars, 'HorizontalAlignment', 'center', 'FontSize', 14);
        end
    end

    xticks(1:num_regions);
    xticklabels(ordered_regions);
    xtickangle(0);

    legend({'Left', 'Right'});
    ylabel('Δz-score');
    title(sprintf('Left vs Right per region (Δz-score, %s)', int_name{t}));

    y_max_overall = max(max(means + errs));
    y_min_overall = min(min(means - errs));
    y_min = floor(y_min_overall / 0.5) * 0.5;
    y_max = ceil(y_max_overall / 0.5) * 0.5;
    ylim([y_min, y_max]);
    yticks(y_min:0.5:y_max);

    grid on;

    % 저장
    fname = sprintf('Bar_Regionwise_LR_ddfmean_%s', int_name{t});
    set(fig,'renderer','Painters');
    print('-depsc','-tiff','-r300','-painters', fullfile(taskDir, [fname '.eps']));
    print(fig, fullfile(taskDir, [fname '.png']), '-dpng', '-r300');
    saveas(fig, fullfile(taskDir, [fname '.fig']));
    close(fig);
end


%% Filtered data

masked_dataF_1x = zeros(n_ints, n_units);
for idx = 1:n_ints
    props = struct();
    props.int_name = intervals{idx}; % interval_all if want 'resting' state too.
    props.atlas_map = atlas_map;
    props.node_positions = node_positions;
    props.colormap_list = colormap_list; 
    props.colormap = 'custom_PiYG';
    props.y_limit = y_limit;
    props.save_path = fullfile(taskDir, 'filt');
    props.sym_option = true;         
    props.unit_size = 0;  
    props.size_list = [0.5, 50; 1.0, 300; 1.3, 700; 1.7, 1200; 1.7, 1700];
    props.trial_type = '1x';

    if ~exist(props.save_path, 'dir')
        mkdir(props.save_path);
    end
    
    % 데이터 및 마스크 설정
    if strcmp(props.int_name, 'resting')
        unit_mask = true(n_units, 1);
        unit_data = resting_mean(:, idx); 
    else
        unit_mask = h_valuesF1(:, idx) == 1;
        unit_data = ddf_mean(idx, :)';
    end

    masked_dataF_1x(idx, :) = atlas_plot(unit_data, unit_mask, props);
end

masked_dataF_2x = zeros(n_ints, n_units);
for idx = 1:n_ints
    props = struct();
    props.int_name = intervals{idx};
    props.atlas_map = atlas_map;
    props.node_positions = node_positions;
    props.colormap_list = colormap_list; 
    props.colormap = 'custom_PiYG';
    props.y_limit = y_limit;
    props.save_path = fullfile(taskDir, 'filt');
    props.sym_option = true;         
    props.unit_size = 0;  
    props.size_list = [0.5, 50; 1.0, 300; 1.3, 700; 1.7, 1200; 1.7, 1700];
    props.trial_type = '2x';

    if ~exist(props.save_path, 'dir')
        mkdir(props.save_path);
    end
    
    if strcmp(props.int_name, 'resting')
        unit_mask = true(n_units, 1);
        unit_data = resting_mean2(:, idx); 
    else
        unit_mask = h_valuesF2(:, idx) == 1;
        unit_data = ddf_mean2(idx, :)';
    end

    masked_dataF_2x(idx, :) = atlas_plot(unit_data, unit_mask, props);
end
save(fullfile(taskDir,'masked_data.mat'),'masked_dataF_1x','masked_dataF_2x');
winopen(taskDir);
return;
%% 
% ΔF/F 차이 계산: mean(2x - 1x)
mean_diff = zeros(n_ints, n_units);
for eidx = 1:n_ints
    for node = 1:n_units
        vals = zeros(num_animals, 1);
        for aidx = 1:num_animals
            vals(aidx) = avg_zs_2x{aidx}(eidx, node) - avg_zs_1x{aidx}(eidx, node);
        end
        mean_diff(eidx, node) = mean(vals, 'omitnan');
    end
end

%% 
df_sym_diff = zeros(size(ddf_mean));
for idx = 1:n_ints
    interval_name = intervals{idx};
    zs_vals = mean_diff(idx, :)';  % [node x 1]

    if sym_option
        zs_temp = reshape(zs_vals, 2, [])';     % [14 x 2]
        zs_sym = mean(zs_temp, 2);              % [14 x 1]
        zs_data = repelem(zs_sym, 2);           % [28 x 1]
    else
        zs_data = zs_vals;
    end
    df_sym_diff(idx, :) = zs_data';
    % 시각화할 좌표
    x_pos = node_positions(:, 1);
    y_pos = node_positions(:, 2);
    node_color = zs_data;

    % node size (절대값 기반)
    abs_vals = abs(node_color);
    node_size = zeros(size(abs_vals));
    node_size(abs_vals < 0.01) = 100;
    node_size(abs_vals >= 0.01 & abs_vals < 0.1) = 300;
    node_size(abs_vals >= 0.1 & abs_vals < 0.2) = 700;
    node_size(abs_vals >= 0.2 & abs_vals < 0.3) = 1200;
    node_size(abs_vals >= 0.3) = 1700;

    % plot
    figure; hold on;
    axis equal; axis off;
    contour(atlas_map, [0.5 0.5], 'k', 'LineWidth', 1.5);
    scatter(x_pos, y_pos, node_size, node_color, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);

    colormap(custom_PiYG);
    caxis(y_limit_ddf);
    colorbar;
    title(sprintf('ΔZ-score (2x - 1x, %s)', interval_name));
    set(gca, 'YDir', 'reverse');

    % save
    if sym_option, suffix = '_sym'; else, suffix = ''; end
    fbase = sprintf('diff_ddf_atlas_nosig%s_%s_%s%s', suffix, task, interval_name, subfix);
    saveas(gcf, fullfile(taskDir, [fbase, '.png']));
    print('-depsc','-tiff','-r300', '-painters', fullfile(taskDir, [fbase, '.eps']));
end
save(fullfile(taskDir, sprintf('PLSR_symdiff_%s%s.mat', task, subfix)), 'df_sym_diff');
%% 
if figure_plot_option == false
    return;
end
%% 

% topFolders에서 마지막 폴더 이름 추출
folderNames = cellfun(@(x) string(x), cellfun(@(x) getLastFolderName(x), topFolders, 'UniformOutput', false));

% GUI를 통한 폴더 선택
[selectedIdx, tf] = listdlg('PromptString', 'Select one folder:', ...
                            'ListString', folderNames, ...
                            'SelectionMode', 'single');

if ~tf
    disp('No folder selected. Exiting...');
    return;
end

% 선택된 폴더 이름과 index
selectedFolder = folderNames{selectedIdx};
fprintf('Selected folder: %s\n', selectedFolder);

% Trial 번호 입력 받기
answer = inputdlg('Enter trial number (0 for random):', 'Trial Selection', [1 50], {'0'});
trial_num = str2double(answer{1});
if isnan(trial_num)
    error('Invalid input. Trial number must be numeric.');
end
topFolder = topFolders{selectedIdx};
subFolders = dir(topFolder);
subFolders = subFolders([subFolders.isdir]);
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));
[~, idx_recent] = max([subFolders.datenum]);
recentFolder = fullfile(topFolder, subFolders(idx_recent).name);
fc_file = fullfile(recentFolder, 'imaging', 'df_data', 'FC_data_all.mat');
if ~isfile(fc_file)
    error('File not found: %s', fc_file);
end

fc_data = load(fc_file);
if ~isfield(fc_data, 'zscore_stacked')
    error('zscored_stacked variable not found in %s', fc_file);
end

zscore_stacked = fc_data.zscore_stacked;

% select trial number
num_trials = size(zscore_stacked, 2);
if trial_num == 0
    rng('shuffle');
    trial_num = randi(num_trials);
    fprintf('Random trial selected: %d\n', trial_num);
elseif trial_num > num_trials || trial_num < 1
    error('Invalid trial number. Must be between 1 and %d.', num_trials);
end

% extract selected trial
trial_data = zscore_stacked{:, trial_num};
fprintf('Loaded trial #%d from %s\n', trial_num, selectedFolder);
%% 
% 옵션 설정
even_only = true;  % ← 여기서 옵션만 바꾸면 됩니다
color_option = true;

% 노드 선택
if even_only
    node_indices = 2:2:size(trial_data, 2);  % 짝수 노드만
else
    node_indices = 1:size(trial_data, 2);    % 전체 노드
end
% 프레임 관련 기본 변수
num_frames = size(trial_data, 1);
onset_frame = 90;
fRate = 20;  % Hz

% x축 시간 (초)로 변환
time_axis = (1:num_frames) - onset_frame;
time_axis = time_axis / fRate;

% 2초 간격으로 xtick 설정 (예: -4, -2, 0, 2, 4, ..., 14)
min_time = floor(min(time_axis) / 2) * 2;
max_time = ceil(max(time_axis) / 2) * 2;
xtick_values = min_time:2:max_time;

% xtick에 해당하는 frame 위치 찾기
xtick_frames = round(xtick_values * fRate + onset_frame);

% 데이터 추출 및 z-score
data = trial_data(:, node_indices);
% data = zscore(data);  % 각 노드별로 정규화

% 노드 수 및 정렬
[num_frames, n_units] = size(data);
gap = 4;  % 노드 간 y 간격
offset = (0:n_units-1) * gap;

% 노드 순서 반전 (작은 번호가 위)
reorder = n_units:-1:1;
data = data(:, reorder);
node_indices_flipped = node_indices(reorder);

if color_option
    all_colors = hsv(size(trial_data, 2));  % 모든 노드용 색상
    line_colors = all_colors(node_indices_flipped, :);  % 현재 사용된 노드에 대응
    xline_color = 'k';  % 검정색
else
    line_colors = repmat([0 0 0], n_units, 1);  % 검정색 반복
    xline_color = 'r';  % 빨간색
end

% 데이터에 offset 적용
shifted_data = data + offset;

% 플롯 시작
figure; hold on;
for i = 1:n_units
    plot(shifted_data(:, i), 'Color', line_colors(i, :), 'LineWidth', 1.2);
end

% motor on/off 선 추가
if color_option
    
else
   
end

% motor on/off 표시
xline(90, '--', 'on', ...
    'Color', xline_color, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'center', ...
    'LineWidth', 1.5);

xline(290, '--', 'off', ...
    'Color', xline_color, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'center', ...
    'LineWidth', 1.5);

% 축 설정
xlabel('Frame');
ylabel('Activity (offset by node)');
title('Node Activity (Z-scored)');
set(gca, 'XTick', xtick_frames, 'XTickLabel', string(xtick_values));
xlabel('Time (s)');
set(gca, 'YTick', offset, ...
         'YTickLabel', strcat("Node ", string(node_indices_flipped)));
xlim([1, num_frames]);
ylim([min(offset)-gap, max(offset)+gap]);
set(gca, 'YDir', 'normal');

% z-score 스케일 바 설정
scale_length = 2.5;  % 예: 2.5 z-score
bar_x = num_frames + 10;  % x 위치 (오른쪽 약간 바깥)
bar_y_start = offset(1);  % 제일 위에 있는 trace 위치 기준

% 플롯에 추가
line([bar_x, bar_x], [bar_y_start, bar_y_start + scale_length], ...
    'Color', 'k', 'LineWidth', 2)

text(bar_x + 5, bar_y_start + scale_length/2, ...
     sprintf('%.1f z-score', scale_length), ...
     'Rotation', 90, 'HorizontalAlignment', 'center')

% ⬅️ xlim, ylim 확대해서 바 안 잘리게
xlim([1, bar_x + 30]);
ylim([min(offset)-gap, max(offset)+gap + scale_length]);
set(gcf, 'windowstate', 'maximized')

saveas(gcf, fullfile(taskDir, sprintf('example_node_plot_%s_%s.png', selectedFolder, num2str(trial_num))));
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', fullfile(taskDir, sprintf('example_node_plot_%s_%s.eps', selectedFolder, num2str(trial_num))));
