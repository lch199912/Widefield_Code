clear; clc; close all;
%% Base parameters
node_selected = 1:2:28; % select odd node idx to plot
y_limit_dR = [-0.4, 0.4]; 
fix_node_size = 00;
plot_option = true;

fisher_transform = true;
paired_option = true; alpha = 0.1;

sym_option = true;
%% Choose task
task_list = {'treadmill', 'wheel', 'disk'};
task_choice = questdlg('Select Task:', 'Task Selection', ...
    'treadmill', 'wheel', 'disk', 'treadmill');
if isempty(task_choice) 
    disp('Task selection canceled. Exiting...');
    return;
end
task_map = containers.Map({'treadmill', 'wheel', 'disk'}, {'treadmill', 'wheel', 'disk'});
task = task_map(task_choice);

disp(['Selected Task: ', task]);

%% Load data
resultDir = fullfile('D:\data analysis\results', 'PLSR_dR_all'); % modify this raw
taskDir = fullfile(resultDir, task);
if ~exist(taskDir, 'dir')
    mkdir(taskDir);
end

file_name = sprintf('PLSR_raw_%s_zscored.mat', task);
mean_data = load(fullfile(taskDir, file_name));
dR_avg1 = mean_data.avg_dR_1x; dR_avg2 = mean_data.avg_dR_2x;
ZdR_avg1 = mean_data.avgZ_dR_1x; ZdR_avg2 = mean_data.avgZ_dR_2x;
ZR_avg1 = mean_data.avgZ_R_1x; ZR_avg2 = mean_data.avgZ_R_2x;

fname = sprintf('PLSR_mean_%s_zscored.mat', task);
mdata = load(fullfile(taskDir, fname));
dR_mean1 = mdata.mean_dR_1x; dR_mean2 = mdata.mean_dR_2x;

load('D:\reference\colormap_data.mat');
load('D:\reference\sorted_node_positions.mat');
load('D:\reference\reference_atlas.mat');
%% 
intervals = {'start', 'continue', 'stop'}; interval_all = [{'resting'}, intervals];
n_ints = length(intervals);
num_animals = size(dR_avg1, 1);
n_units = size(dR_avg1{1}, 1);

dR_interval1 = cell(n_ints, 1); dR_interval2 = cell(n_ints, 1);
for eidx = 1:n_ints
    dR_interval1{eidx} = NaN(n_units, n_units, num_animals); dR_interval2{eidx} = NaN(n_units, n_units, num_animals);
    for iidx = 1:num_animals
        dR_interval1{eidx}(:,:,iidx) = dR_avg1{iidx}(:,:,eidx); dR_interval2{eidx}(:,:,iidx) = dR_avg2{iidx}(:,:,eidx);
    end
end

%% t-test validation
% Initialize results
p_values1 = NaN(n_units, n_units, n_ints);
h_values1 = NaN(n_units, n_units, n_ints);
all_p_values1 = []; 
all_indices1 = []; 

% Perform one-sample t-test on ΔZ (across animals)
for eidx = 1:n_ints
    for row = 1:n_units
        for col = 1:n_units
            deltaZ1 = zeros(num_animals, 1); % ΔZ = Z_event - Z_resting
            
            for aidx = 1:num_animals
                deltaZ1(aidx) = ZdR_avg1{aidx}(row, col, eidx); 
            end

            % 자기자신 (diagonal)은 제외
            if row == col
                continue;
            end

            if paired_option
                % Paired-sample t-test
                eventZ = zeros(num_animals, 1);
                restZ  = zeros(num_animals, 1);
                for aidx = 1:num_animals
                    eventZ(aidx) = ZR_avg1{aidx}(row, col, eidx + 1);
                    restZ(aidx)  = ZR_avg1{aidx}(row, col, 1);
                end
                [~, p] = ttest(eventZ, restZ, 'Alpha', alpha);

            else
                % One-sample t-test on ΔZ
                deltaZ = zeros(num_animals, 1);
                for aidx = 1:num_animals
                    deltaZ(aidx) = ZdR_avg1{aidx}(row, col, eidx);  % ΔZ = Z_event - Z_resting
                end
                [~, p] = ttest(deltaZ, 0, 'Alpha', alpha);
            end

            p_values1(row, col, eidx) = p;
            all_p_values1 = [all_p_values1; p];
            all_indices1 = [all_indices1; row, col, eidx];
        end
    end
end

% Initialize results
p_values2 = NaN(n_units, n_units, n_ints);
h_values2 = NaN(n_units, n_units, n_ints);
all_p_values2 = []; 
all_indices2 = []; 

% Perform one-sample t-test on ΔZ (across animals)
for eidx = 1:n_ints
    for row = 1:n_units
        for col = 1:n_units
            deltaZ2 = zeros(num_animals, 1); % ΔZ = Z_event - Z_resting
            
            for aidx = 1:num_animals
                deltaZ2(aidx) = ZdR_avg2{aidx}(row, col, eidx); 
            end

            % 자기자신 (diagonal)은 제외
            if row == col
                continue;
            end
            if paired_option
                % Paired-sample t-test
                eventZ = zeros(num_animals, 1);
                restZ  = zeros(num_animals, 1);
                for aidx = 1:num_animals
                    eventZ(aidx) = ZR_avg1{aidx}(row, col, eidx + 1);
                    restZ(aidx)  = ZR_avg1{aidx}(row, col, 1);
                end
                [~, p] = ttest(eventZ, restZ, 'Alpha', alpha);

            else
                % One-sample t-test on ΔZ
                deltaZ = zeros(num_animals, 1);
                for aidx = 1:num_animals
                    deltaZ(aidx) = ZdR_avg1{aidx}(row, col, eidx);  % ΔZ = Z_event - Z_resting
                end
                [~, p] = ttest(deltaZ, 0, 'Alpha', alpha);
            end

            p_values2(row, col, eidx) = p;
            all_p_values2 = [all_p_values2; p];
            all_indices2 = [all_indices2; row, col, eidx];
        end
    end
end


%% Apply FDR correction
[h_fdr1, crit_p1] = fdr_bh(all_p_values1, alpha); % FDR correction 적용

% 유의한 index 정리
significant_idx1 = find(h_fdr1 == 1);
sorted_idx1 = all_indices1(significant_idx1, :); % 유의한 row, col, interval index 정렬
sorted_p_values1 = all_p_values1(significant_idx1); % 유의한 p-value 정렬

% h_values에 저장
for i = 1:length(significant_idx1)
    row = sorted_idx1(i, 1);
    col = sorted_idx1(i, 2);
    eidx = sorted_idx1(i, 3);
    h_values1(row, col, eidx) = 1;
end

[h_fdr2, crit_p2] = fdr_bh(all_p_values2, alpha); % FDR correction 적용

% 유의한 index 정리
significant_idx2 = find(h_fdr2 == 1);
sorted_idx2 = all_indices2(significant_idx2, :); % 유의한 row, col, interval index 정렬
sorted_p_values2 = all_p_values2(significant_idx2); % 유의한 p-value 정렬

% h_values에 저장
for i = 1:length(significant_idx2)
    row = sorted_idx2(i, 1);
    col = sorted_idx2(i, 2);
    eidx = sorted_idx2(i, 3);
    h_values2(row, col, eidx) = 1;
end
%% save t-test data
parameters = struct();
parameters.node_selected = node_selected;
parameters.fisher_transform = fisher_transform;
parameters.plot_option = plot_option;
parameters.num_intervals = n_ints;
parameters.num_animals = num_animals;
parameters.num_nodes = n_units;


%% diff t-test
% 옵션 설정
use_baseline = true;  % true: (Z2 - Zrest2) - (Z1 - Zrest1), false: Z2 - Z1

% 초기화
p_values_diff = NaN(n_units, n_units, n_ints);
h_values_diff = NaN(n_units, n_units, n_ints);
all_p_values_diff = [];
all_indices_diff = [];

% 조건별 t-test
for eidx = 1:n_ints
    for row = 1:n_units
        for col = 1:n_units
            if row == col
                continue;  % 자기 자신은 제외
            end

            x1 = zeros(num_animals, 1);  % 1x 조건
            x2 = zeros(num_animals, 1);  % 2x 조건

            for aidx = 1:num_animals
                if use_baseline
                    % ΔZ = Z_event - Z_resting
                    z1 = ZdR_avg1{aidx}(row, col, eidx);
                    z2 = ZdR_avg2{aidx}(row, col, eidx);
                else
                    % 원본 Z값에서 비교 (Z_event 값 자체)
                    z1 = Z_event1{aidx}(row, col, eidx);  % ← 미리 만들어져 있어야 함
                    z2 = Z_event2{aidx}(row, col, eidx);  % ← 미리 만들어져 있어야 함
                end
                x1(aidx) = z1;
                x2(aidx) = z2;
            end

            % paired t-test
            [~, p] = ttest(x2, x1, 'Alpha', alpha);
            p_values_diff(row, col, eidx) = p;

            all_p_values_diff = [all_p_values_diff; p];
            all_indices_diff = [all_indices_diff; row, col, eidx];
        end
    end
end

% FDR 보정
[h_fdr_diff, crit_p_diff] = fdr_bh(all_p_values_diff, alpha);

% 유의성 저장
for i = 1:length(h_fdr_diff)
    if h_fdr_diff(i)
        row = all_indices_diff(i, 1);
        col = all_indices_diff(i, 2);
        eidx = all_indices_diff(i, 3);
        h_values_diff(row, col, eidx) = 1;
    end
end
%% 
save(fullfile(taskDir, sprintf('ttest_atlas_%s.mat', task)), 'sorted_idx1', 'sorted_p_values1', 'h_values1', 'p_values1', 'crit_p1', ...
    'sorted_idx2', 'sorted_p_values2', 'h_values2', 'p_values2', 'crit_p2','parameters', 'h_values_diff');
disp('t-test done... saved t-test data...');
%% skip plot
if plot_option == false
    disp('plot_option is false. skip plotting network graph.');
    return;
end


x_pos = node_positions(:, 1);
y_pos = node_positions(:, 2);
% Define left (odd) and right (even) hemisphere nodes
left_nodes = 1:2:n_units; 
right_nodes = 2:2:n_units; 
%% average contra ipsi
masked_1x = zeros(size(dR_mean1)); masked_2x = zeros(size(dR_mean2));
dR_mean_sym1 = NaN(n_ints, n_units);
dR_mean_sym2 = NaN(n_ints, n_units);
mean_diff_sym = NaN(n_ints, n_units);
dR_masked_1x = NaN(n_ints, n_units);
dR_masked_2x = NaN(n_ints, n_units);
mean_diff_sym_raw = NaN(n_ints, n_units);
dR_unmasked_1x = NaN(n_ints, n_units);
dR_unmasked_2x = NaN(n_ints, n_units);

for eidx = 1:n_ints
    h_mask = h_values1(:,:,eidx);
    curr1 = dR_mean1(:,:,eidx);
    h_mask2 = h_values2(:,:,eidx);
    curr2 = dR_mean2(:,:,eidx);
    diff = (curr2 - curr1) ./ curr1;

    % ✅ 비마스킹 평균 (대각 제외)
    curr1_no_diag = curr1;
    curr1_no_diag(logical(eye(n_units))) = NaN;
    dR_unmasked_1x(eidx, :) = mean(curr1_no_diag, 2, 'omitnan')';

    curr2_no_diag = curr2;
    curr2_no_diag(logical(eye(n_units))) = NaN;
    dR_unmasked_2x(eidx, :) = mean(curr2_no_diag, 2, 'omitnan')';

    % ✅ 마스킹된 ΔR
    curr1(logical(eye(n_units))) = NaN;
    curr1(isnan(h_mask)) = NaN;
    masked_1x(:,:, eidx) = curr1;
    dR_masked_1x(eidx, :) = mean(curr1, 2, 'omitnan')';

    curr2(logical(eye(n_units))) = NaN;
    curr2(isnan(h_mask2)) = NaN;
    masked_2x(:,:, eidx) = curr2;
    dR_masked_2x(eidx, :) = mean(curr2, 2, 'omitnan')';

    % ✅ 좌우 평균 및 차이
    for i = 1:length(left_nodes)
        li = left_nodes(i);
        ri = right_nodes(i);
        sym1_val = mean([curr1(li, :), curr1(ri, :)], 'all', 'omitnan');
        sym2_val = mean([curr2(li, :), curr2(ri, :)], 'all', 'omitnan');
        diff_val = sym2_val - sym1_val;
        diff_sym = mean([diff(li, :), diff(ri, :)], 'all','omitnan');
        dR_mean_sym1(eidx, [li ri]) = sym1_val;
        dR_mean_sym2(eidx, [li ri]) = sym2_val;
        mean_diff_sym(eidx, [li ri]) = diff_val;
        mean_diff_sym_raw(eidx, [li ri]) = diff_sym;
    end
end
save(fullfile(taskDir, sprintf('PLSR_symFC_%s.mat', task)), 'dR_mean_sym1', 'dR_mean_sym2','dR_masked_1x','dR_masked_2x','mean_diff_sym','mean_diff_sym_raw');
% return;
%% 
% node_size(abs(dR_val) < 0.01) = 100;     % 작은 값 (0~0.01)
% node_size(abs(dR_val) >= 0.01 & abs(dR_val) < 0.03) = 300;  % 중간 값 (0.01~0.03)
% node_size(abs(dR_val) >= 0.03 & abs(dR_val) < 0.05) = 700;  % 큰 값 (0.03~0.05)
% node_size(abs(dR_val) >= 0.05) = 1200;   % 매우 큰 값 (>=0.05)

%% 1x delta R
for idx = 1:n_ints
    interval_name = intervals{idx};

    % 데이터 선택 및 전처리
    if sym_option
        dR_data = dR_mean_sym1(idx, :)';
    else
        dR_data = dR_masked_1x(idx, :)';
    end
    dR_data(isnan(dR_data)) = 0;
    unit_mask = dR_data ~= 0;

    % properties 구조체 정의
    properties = struct();
    properties.sym_option     = sym_option;
    properties.fix_node_size  = exist('fix_node_size', 'var') && fix_node_size ~= 0;
    properties.unit_size      = 0;
    properties.node_positions = node_positions;
    properties.atlas_map      = atlas_map;
    properties.save_path      = taskDir;
    properties.int_name       = interval_name;
    properties.trial_type     = sprintf('dR_atlas_1x_%s', task);
    properties.y_limit        = y_limit_dR;
    properties.colormap_list  = struct('RdBu', custom_RdBu);  % ✅ custom colormap 구조체화
    properties.colormap       = 'RdBu';  % 사용할 필드명
    properties.size_list = [          % size_list는 abs ΔR 기준
        0.05, 50;
        0.15, 200;
        0.3,  600;
        0.5, 1000;
        1.0, 1500
    ];

    % 함수 실행
    atlas_plot(dR_data, unit_mask, properties);
end
%% 

% 영역별 유닛 인덱스 정의 (1~28 유닛 기준)
region_indices = {
    1:6;     % M2
    7:10;    % M1
    11:16;   % S1
    17:18;   % Aud
    19:24;   % Vis
    25:28    % RSC
};

n_intervals = size(dR_unmasked_1x, 1);
n_regions = numel(region_indices);
dR_region1 = NaN(n_intervals, n_regions);  % 결과 저장용 (3×6)

% 평균 계산
for i = 1:n_intervals
    for r = 1:n_regions
        unit_idx = region_indices{r};
        region_values = dR_unmasked_1x(i, unit_idx);
        dR_region1(i, r) = mean(region_values, 'omitnan');  % NaN 무시
    end
end
%% 
num_animals = numel(ZdR_avg1);
n_units = 28;
n_intervals = 3;
ZdR_selfavg1 = cell(num_animals, 1);

for i = 1:num_animals
    mat = ZdR_avg1{i};  % [28×28×3]
    
    % 열 방향 평균 (각 unit → 다른 unit들과의 평균 연결성)
    avg_over_columns = mean(mat, 2, 'omitnan');  % [28×1×3]
    
    ZdR_selfavg1{i} = squeeze(avg_over_columns);  % → [28×3]
end
% 영역별 유닛 인덱스 정의
region_indices = {
    1:6;     % M2
    7:10;    % M1
    11:16;   % S1
    17:18;   % Aud
    19:24;   % Vis
    25:28    % RSC
};

n_animals = numel(ZdR_selfavg1);
n_intervals = size(ZdR_selfavg1{1}, 2);  % = 3
n_regions = numel(region_indices);

ZdR_region1 = cell(n_animals, 1);  % 결과 저장용: 5x1 cell, 각 cell은 [3 x 6]

for i = 1:n_animals
    data = ZdR_selfavg1{i};  % [28 x 3]
    region_avg = NaN(n_intervals, n_regions);  % [3 x 6]

    for r = 1:n_regions
        unit_idx = region_indices{r};
        region_values = data(unit_idx, :);  % [n_units_in_region x 3]
        region_avg(:, r) = mean(region_values, 1, 'omitnan');  % [1 x 3] → [3 x 1] transpose
    end

    ZdR_region1{i} = region_avg;  % [3 x 6]
end

save(fullfile(taskDir, ['dR_region_mean_', task, '.mat']), 'ZdR_selfavg1', 'ZdR_region1', 'dR_region1', 'dR_unmasked_1x');
winopen(taskDir); return;
%% 2x delta R
% for idx = 1:num_intervals
%     interval_name = intervals{idx};
% 
%     if sym_option == true
%         dR_data2 = dR_mean_sym2(idx, :);
%     else
%         dR_data2 = dR_masked_2x(idx, :);
%     end
%     dR_data2(isnan(dR_data2)) = 0;
%     sig_idx = dR_data2 ~= 0;
%     x_pos = node_positions(sig_idx, 1);
%     y_pos = node_positions(sig_idx, 2);
%     node_color = dR_data2(sig_idx);
% 
%     % ✅ 구간별 node size 설정
%     abs_vals = abs(node_color);
%     node_size = zeros(size(abs_vals));
%     node_size(abs_vals < 0.1) = 100;
%     node_size(abs_vals >= 0.1 & abs_vals < 0.2) = 200;
%     node_size(abs_vals >= 0.2 & abs_vals < 0.3) = 600;
%     node_size(abs_vals >= 0.3 & abs_vals < 0.45) = 1000;
%     node_size(abs_vals >= 0.45) = 1500;
% 
%     % 시각화
%     figure;    
%     set(gcf, 'Color', 'w');
%     axis equal;
%     axis off;
%     hold on;
%     contour(atlas_map, [0.5 0.5], 'k', 'LineWidth', 1.5);
%     scatter(x_pos, y_pos, node_size, node_color, 'filled', ...
%         'MarkerEdgeColor', 'k', 'LineWidth', 1);
% 
% %     colormap(flipud(slanCM('spectral')));
%     colormap(custom_RdBu);
% %     colormap('jet'); 
%     caxis(y_limit_dR);
%     colorbar;
%     title(sprintf('mean ΔR (2x, %s)', interval_name));
%     set(gca, 'YDir', 'reverse'); 
%     hold off;
% 
%     % save figure
%     if sym_option
%         suffix = '_sym';
%     else
%         suffix = '';
%     end
%     if exist('fix_node_size', 'var') && fix_node_size ~= 0
%         plot_fname = fullfile(taskDir, sprintf('fixed_dR_atlas_2x%s_%s_%s.png', suffix, task, interval_name));
%         plot_epsname = fullfile(taskDir, sprintf('fixed_dR_atlas_2x%s_%s_%s.eps', suffix, task, interval_name));
%     else
%         plot_fname = fullfile(taskDir, sprintf('dR_atlas_2x%s_%s_%s.png', suffix, task, interval_name));
%         plot_epsname = fullfile(taskDir, sprintf('dR_atlas_2x%s_%s_%s.eps', suffix, task, interval_name));
%     end
%     
%     saveas(gcf, plot_fname);
%     set(gcf,'renderer','Painters');
%     print('-depsc','-tiff','-r300', '-painters', plot_epsname);
% end
%% 2x-1x diff delta R
% for idx = 1:num_intervals
%     dR_val = mean_diff_sym_raw(idx, :);  % [1 x num_nodes], 좌우 평균 ΔR (2x - 1x)
% 
%     % 🔍 FC 변화가 유의한 node만 남기기 위해 mask 생성
%     % h_values_diff는 [node x node x interval]
%     h_mask = squeeze(h_values_diff(:, :, idx));     % 현재 interval의 FC 유의성 [node x node]
%     sig_nodes = any(h_mask == 1, 2);                % 이 노드가 유의한 연결을 하나라도 가지면 true
% 
%     % 🔍 마스킹된 ΔR 값 만들기
%     dR_filtered = dR_val;
%     dR_filtered(~sig_nodes') = 0;   % 유의하지 않은 node는 0으로 masking
% 
%     % 시각화 시작
%     figure;
%     set(gcf, 'Color', 'w'); axis equal; axis off; hold on;
% 
%     contour(atlas_map, [0.5 0.5], 'k', 'LineWidth', 1.5);
% 
%     node_size = abs(dR_filtered) * 10000;
%     node_color = dR_filtered;
% 
%     scatter(x_pos, y_pos, node_size, node_color, 'filled', ...
%         'MarkerEdgeColor', 'k', 'LineWidth', 1); 
% 
%     for i = 1:num_nodes
%         text(x_pos(i), y_pos(i), num2str(i), 'Color', 'w', ...
%              'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
%     end
% 
%     colormap(custom_RdBu);
%     caxis([-0.1, 0.1]);
%     colorbar;
%     title(sprintf('ΔR Avg 2x-1x - %s (Significant Only)', intervals{idx}));
%     set(gca, 'YDir', 'reverse');
% 
%     fname = fullfile(taskDir, sprintf('dR_sym_nodewise_diff_filtered_%s.png', intervals{idx}));
%     epsname = fullfile(taskDir, sprintf('dR_sym_nodewise_diff_filtered_%s.eps', intervals{idx}));
%     saveas(gcf, fname);
%     print('-depsc','-tiff','-r300', '-painters', epsname);
% end
