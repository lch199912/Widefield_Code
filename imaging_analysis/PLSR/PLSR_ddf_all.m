clear;clc;close all;
%% parameters
resting_time = [11, 70]; % resting time index with second (0 = motor onset)
walking_time = [91, 150; 191, 250; 291, 350]; % time to compare (0 = motor onset)

calculate_type = 'zscore';
%% Load ICA file for df extract
addpath(genpath('E:\Final Codes'))
analysis_dir = 'D:\data analysis\locomotion training\'; cd(analysis_dir);
m_analysis_dir = uigetdir(analysis_dir, 'Select Mouse Folder');
if m_analysis_dir == 0
    error('No directory selected. Exiting...');
end
[analysis_base_dir, mouse] = fileparts(m_analysis_dir);
fprintf('Selected Mouse: %s\n', mouse);

day_folders = dir(fullfile(m_analysis_dir)); 
day_folders = day_folders([day_folders.isdir]);
day_folders = day_folders(~ismember({day_folders.name}, {'.', '..'}));
if isempty(day_folders)
    error('No day folders found in the selected mouse folder.');
end
fprintf('Available Days for Mouse %s:\n', mouse);
for i = 1:length(day_folders)
    fprintf('%d: %s\n', i, day_folders(i).name);
end
day = day_folders(1).name; 
fprintf('Selected Day: %s\n', day);

fdir = fullfile(analysis_dir, mouse, day, 'imaging');
load(fullfile(fdir, "atlas.mat"));
brainmask = ~isnan(brainmask); % convert into boolean
load(fullfile(fdir, 'locaNMF','properties.mat'));
dfDir = fullfile(fdir, 'df_data', 'PLSR_all');
if strcmp(calculate_type, 'raw')
    load(fullfile(dfDir, 'PLSR_z-score.mat'), 'filtered_data')
    subfix = '';
elseif strcmp(calculate_type, 'zscore')
    load(fullfile(dfDir, 'PLSR_z-score_zscored.mat'), 'filtered_data')
    subfix = '_zscored';
else
    error('wrong type')
end
%% 

zscore_all = filtered_data;
zscore_1x_all = filtered_data(list_1x); zscore_2x_all = filtered_data(list_2x); 
zscore_1x = mean_cell_matrices(zscore_1x_all);
zscore_2x = mean_cell_matrices(zscore_2x_all);

%% 

num_nodes = size(zscore_1x, 2);
num_intervals = size(walking_time, 1);

zs_rest_1x = mean(zscore_1x(resting_time(1):resting_time(2), :), 1, 'omitnan'); 
zs_rest_2x = mean(zscore_2x(resting_time(1):resting_time(2), :), 1, 'omitnan'); 

mean_dzs_1x = zeros(num_intervals, num_nodes);
mean_dzs_2x = zeros(num_intervals, num_nodes);
mean_zs_1x = zeros(num_intervals, num_nodes);
mean_zs_2x = zeros(num_intervals, num_nodes);

for iidx = 1:num_intervals
    zs_mean_1x = mean(zscore_1x(walking_time(iidx, 1):walking_time(iidx, 2), :), 1, 'omitnan');
    zs_mean_2x = mean(zscore_2x(walking_time(iidx, 1):walking_time(iidx, 2), :), 1, 'omitnan');

    mean_zs_1x(iidx, :) = zs_mean_1x;
    mean_zs_2x(iidx, :) = zs_mean_2x;
    mean_dzs_1x(iidx, :) = zs_mean_1x - zs_rest_1x;
    mean_dzs_2x(iidx, :) = zs_mean_2x - zs_rest_2x;
end
%% 
% save Mean Δ df 
% save(fullfile(dfDir, 'PLSR_ddf.mat'), 'zs_rest_1x', 'zs_rest_2x', 'mean_zs_1x', 'mean_zs_2x', 'mean_dzs_1x', 'mean_dzs_2x', 'calculate_type','resting_time','walking_time');
save(fullfile(dfDir, ['PLSR_ddf', subfix, '.mat']), 'zscore_1x', 'zscore_2x','zscore_2x_all', 'zscore_1x_all','zs_rest_1x', 'zs_rest_2x', 'mean_zs_1x', 'mean_zs_2x', 'mean_dzs_1x', 'mean_dzs_2x', 'calculate_type','resting_time','walking_time');
disp('Mean Delta df calculated and saved.');
return;
%% 
% node_labels = 1:28;
% plot_activity_heatmap(zscore_1x, 20, node_labels, 90, 290, [-3, 4], 'jet', 'zscored Heatmap for 1x trial');
% set(gcf, 'WindowState','maximized')
% saveas(gcf, fullfile(dfDir, ['zscore_heatmap_1x', subfix, '.png'])); 
% close(gcf);
% 
% plot_activity_heatmap(zscore_2x, 20, node_labels, 90, 290, [-3, 4], 'jet', 'zscored Heatmap for 2x trial');
% set(gcf, 'WindowState','maximized')
% saveas(gcf, fullfile(dfDir, ['zscore_heatmap_2x', subfix, '.png'])); 
% close(gcf);
% disp('df heatmap saved.')
%% 
data = mean_dzs_1x;
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
    xlim(axis_limits(t, 1:2));
    ylim(axis_limits(t, 3:4));

    xline = linspace(axis_limits(t, 1), axis_limits(t, 2), 100);
    plot(xline, xline, 'k--', 'LineWidth', 1.2);  % 검정 점선

    % 회귀 분석
    mdl = fitlm(x', y');  % row vector -> column vector
    
    % 회귀선 시각화
    x_fit = linspace(axis_limits(t,1), axis_limits(t,2), 100);
    y_fit = predict(mdl, x_fit');
    plot(x_fit, y_fit, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);  % 회색 실선
    
    % 상관계수 계산
    R = corrcoef(x, y);  % 2x2 행렬
    r_value = R(1, 2);
    
    % R² 표시
    text_x = axis_limits(t,1) + 0.05;
    text_y = axis_limits(t,4) - 0.05;
    text(text_x, text_y, sprintf('R² = %.3f', mdl.Rsquared.Ordinary), ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    % (옵션) r-value도 함께 표시
    text(text_x, text_y - 0.05, sprintf('r = %.3f', r_value), ...
        'FontSize', 10, 'FontWeight', 'normal');


    grid on;

    % 저장
    filename_base = sprintf('Scatter_LeftRight_%s', int_name{t});
    saveas(fig, fullfile(dfDir, [filename_base, '.png']));
    saveas(fig, fullfile(dfDir, [filename_base, '.fig']));
    close(fig);
end



