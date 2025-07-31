clear; clc; close all;
%% Base parameters
y_limit_df = [-2.0, 2.0];
fix_node_size = 00; % fix atlas plot node size. if not to use, set this value 0. recommend 500
sym_option = true;
alpha = 0.05;

intervals = {'initiation', 'execution', 'completion'}; interval_all = [{'resting'}, intervals];
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

resultDir = fullfile('D:\data analysis\results', 'avg_ddf');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
taskDir = fullfile(resultDir, task);
if ~exist(taskDir, 'dir')
    mkdir(taskDir);
end
load('D:\reference\colormap_data.mat');
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

for i = 1:numFolders
    topFolder = topFolders{i};
    [~, mouseName] = fileparts(topFolder);
    subFolders = dir(topFolder);
    subFolders = subFolders([subFolders.isdir]);
    subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

    [~, idx_recent] = max([subFolders.datenum]);
    recentFolder = fullfile(topFolder, subFolders(idx_recent).name);

    targetFile = fullfile(recentFolder, 'imaging', 'df_data', 'FC_ddf.mat');

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
        else
            warning('data not found in file %s', targetFile);
        end
    else
        warning('File does not exist: %s', targetFile);
    end
end

%% Calculate mean and std for each dataset
ddf_mean = mean_cell_matrices(avg_dzs_1x); ddf_mean2 = mean_cell_matrices(avg_dzs_2x);
df_mean = mean_cell_matrices(avg_zs_1x); df_mean2 = mean_cell_matrices(avg_zs_2x);
resting_mean = mean_cell_matrices(resting_1x); resting_mean2 = mean_cell_matrices(resting_2x);
ddf_std = std_cell_matrices(avg_dzs_1x); ddf_se = sem_cell_matrices(avg_dzs_1x); 

load('D:\reference\sorted_node_positions.mat');
load('D:\reference\reference_atlas.mat');
save(fullfile(taskDir, 'raw_data.mat'),'avg_dzs_1x','avg_dzs_2x', 'avg_zs_1x','avg_zs_2x', 'resting_1x','resting_2x', 'ddf_mean','ddf_mean2', ...
    'df_mean','df_mean2',"resting_mean",'resting_mean2','ddf_std','ddf_se');return;
%% ΔF/F 값에 대한 t-test 수행 (resting 기준)
% ΔF = df_mean - resting_mean (animal 단위)

% 가정: avg_zs_1x = {N x 1 cell} → 각 cell: [interval x nodes]
%       resting_1x = {N x 1 cell} → 각 cell: [1 x nodes]
% 기본 설정
num_animals = numel(avg_zs_1x);
num_intervals = size(avg_zs_1x{1}, 1);
num_nodes = size(avg_zs_1x{1}, 2);

%% 1x trial t-test

p_values_zs1 = NaN(num_nodes, num_intervals);
h_values_zs1 = NaN(num_nodes, num_intervals);
all_p_values_zs1 = [];
all_indices_zs1 = [];

% ΔF 계산 및 t-test
for eidx = 1:num_intervals
    for node = 1:num_nodes
        deltazs = zeros(num_animals, 1);

        for aidx = 1:num_animals
            deltazs(aidx) = avg_zs_1x{aidx}(eidx, node) - resting_1x{aidx}(node);
        end

        [~, p] = ttest(deltazs, 0, 'Alpha', alpha);
        p_values_zs1(node, eidx) = p;

        all_p_values_zs1 = [all_p_values_zs1; p];
        all_indices_zs1 = [all_indices_zs1; node, eidx];
    end
end

% FDR 보정 적용
[h_fdr_zs1, crit_p_zs1] = fdr_bh(all_p_values_zs1, alpha); 

% 유의한 index 정리
significant_idx_zs1 = find(h_fdr_zs1 == 1);
sorted_idx_zs1 = all_indices_zs1(significant_idx_zs1, :);
sorted_p_values_zs1 = all_p_values_zs1(significant_idx_zs1);

% 유의성 표시 저장
for i = 1:length(significant_idx_zs1)
    node = sorted_idx_zs1(i, 1);
    eidx = sorted_idx_zs1(i, 2);
    h_values_zs1(node, eidx) = 1;
end
%% 2x trial t-test

p_values_zs2 = NaN(num_nodes, num_intervals);
h_values_zs2 = NaN(num_nodes, num_intervals);
all_p_values_zs2 = [];
all_indices_zs2 = [];

% ΔF 계산 및 t-test
for eidx = 1:num_intervals
    for node = 1:num_nodes
        deltazs = zeros(num_animals, 1);

        for aidx = 1:num_animals
            deltazs(aidx) = avg_zs_2x{aidx}(eidx, node) - resting_2x{aidx}(node);
        end

        [~, p] = ttest(deltazs, 0, 'Alpha', alpha);
        p_values_zs2(node, eidx) = p;

        all_p_values_zs2 = [all_p_values_zs2; p];
        all_indices_zs2 = [all_indices_zs2; node, eidx];
    end
end

% FDR 보정 적용
[h_fdr_zs2, crit_p_zs2] = fdr_bh(all_p_values_zs2, alpha); 

% 유의한 index 정리
significant_idx_zs2 = find(h_fdr_zs2 == 1);
sorted_idx_zs2 = all_indices_zs2(significant_idx_zs2, :);
sorted_p_values_zs2 = all_p_values_zs2(significant_idx_zs2);

% 유의성 표시 저장
for i = 1:length(significant_idx_zs2)
    node = sorted_idx_zs2(i, 1);
    eidx = sorted_idx_zs2(i, 2);
    h_values_zs2(node, eidx) = 1;
end
%% save data
save(fullfile(taskDir, sprintf('ttest_%s.mat', task)), ...
    'p_values_zs1', 'h_values_zs1', 'sorted_idx_zs1', 'sorted_p_values_zs1', 'crit_p_zs1', ...
    'p_values_zs2', 'h_values_zs2', 'sorted_idx_zs2', 'sorted_p_values_zs2', 'crit_p_zs2');


%% 🎯 h_masked ΔF/F(zscore) atlas data (1x)
df_sym_masked = zeros(size(ddf_mean));
for idx = 1:num_intervals
    interval_name = intervals{idx};

    % mask 선택
    if strcmp(interval_name, 'resting')
        node_mask = true(1, size(resting_mean, 2));  % 모든 노드 사용
        zs_vals = resting_mean;
    else
        node_mask = h_values_zs1(:, idx) == 1;  % 유의한 노드만
        zs_vals = ddf_mean(idx, :)';           % ΔF/F 값
        zs_vals(~node_mask) = 0;
        zs_data = zs_vals;
        if sym_option == true
            zs_temp = reshape(zs_vals, 2, [])';
            zs_sym = mean(zs_temp, 2);
            zs_data =repelem(zs_sym, 2);
        end
    end
    df_sym_masked(idx, :) = zs_data';
    sig_idx = zs_data ~= 0;
    x_pos = node_positions(sig_idx, 1);
    y_pos = node_positions(sig_idx, 2);
    node_color = zs_data(sig_idx);

    % ✅ 구간별 node size 설정
    abs_vals = abs(node_color);
    node_size = zeros(size(abs_vals));
    node_size(abs_vals < 0.5) = 100;
    node_size(abs_vals >= 0.5 & abs_vals < 1) = 300;
    node_size(abs_vals >= 1 & abs_vals < 1.5) = 700;
    node_size(abs_vals >= 1.5 & abs_vals < 2.0) = 1200;
    node_size(abs_vals >= 2.0) = 1700;

    % 시각화
    figure;    
    set(gcf, 'Color', 'w');
    axis equal;
    axis off;
    hold on;
    contour(atlas_map, [0.5 0.5], 'k', 'LineWidth', 1.5);
    scatter(x_pos, y_pos, node_size, node_color, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);

%     colormap(flipud(slanCM('spectral')));
    colormap(custom_PiYG); 
%     colormap('jet'); 
    caxis(y_limit_df);
    colorbar;
    title(sprintf('mean z-score (1x, %s)', interval_name));
    set(gca, 'YDir', 'reverse'); 
    hold off;

    % save figure
    if sym_option
        suffix = '_sym';
    else
        suffix = '';
    end
    if exist('fix_node_size', 'var') && fix_node_size ~= 0
        plot_fname = fullfile(taskDir, sprintf('fixed_ddf_atlas_1x%s_%s_%s.png', suffix, task, interval_name));
        plot_epsname = fullfile(taskDir, sprintf('fixed_ddf_atlas_1x%s_%s_%s.eps', suffix, task, interval_name));
    else
        plot_fname = fullfile(taskDir, sprintf('ddf_atlas_1x%s_%s_%s.png', suffix, task, interval_name));
        plot_epsname = fullfile(taskDir, sprintf('ddf_atlas_1x%s_%s_%s.eps', suffix, task, interval_name));
    end
    
    saveas(gcf, plot_fname);
    set(gcf,'renderer','Painters');
    print('-depsc','-tiff','-r300', '-painters', plot_epsname);
end
save(fullfile(taskDir, sprintf('PLSR_symdf_%s%s.mat', task, suffix)), 'df_sym_masked'); return;
%% 🎯 h_masked ΔF/F(zscore) atlas data (2x)

for idx = 1:num_intervals
    interval_name = intervals{idx};

    % mask 선택
    if strcmp(interval_name, 'resting')
        node_mask = true(1, size(resting_mean, 2));  % 모든 노드 사용
        zs_vals = resting_mean;
    else
        node_mask = h_values_zs2(:, idx) == 1;  % 유의한 노드만
        zs_vals = ddf_mean2(idx, :)';           % ΔF/F 값
        zs_vals(~node_mask) = 0;
        zs_data = zs_vals;
        if sym_option == true
            zs_temp = reshape(zs_vals, 2, [])';
            zs_sym = mean(zs_temp, 2);
            zs_data =repelem(zs_sym, 2);
        end
    end

    sig_idx = zs_data ~= 0;
    x_pos = node_positions(sig_idx, 1);
    y_pos = node_positions(sig_idx, 2);
    node_color = zs_data(sig_idx);

    % ✅ 구간별 node size 설정
    abs_vals = abs(node_color);
    node_size = zeros(size(abs_vals));
    node_size(abs_vals < 0.5) = 100;
    node_size(abs_vals >= 0.5 & abs_vals < 1) = 300;
    node_size(abs_vals >= 1 & abs_vals < 1.3) = 700;
    node_size(abs_vals >= 1.3 & abs_vals < 1.7) = 1200;
    node_size(abs_vals >= 1.7) = 1700;

    % 시각화
    figure;    
    set(gcf, 'Color', 'w');
    axis equal;
    axis off;
    hold on;
    contour(atlas_map, [0.5 0.5], 'k', 'LineWidth', 1.5);
    scatter(x_pos, y_pos, node_size, node_color, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);

%     colormap(flipud(slanCM('spectral')));
    colormap(custom_PiYG); 
%     colormap('jet'); 
    caxis(y_limit_df);
    colorbar;
    title(sprintf('mean z-score (2x, %s)', interval_name));
    set(gca, 'YDir', 'reverse'); 
    hold off;

    % save figure
    if sym_option
        suffix = '_sym';
    else
        suffix = '';
    end
    if exist('fix_node_size', 'var') && fix_node_size ~= 0
        plot_fname2 = fullfile(taskDir, sprintf('fixed_ddf_atlas_2x%s_%s_%s.png', suffix, task, interval_name));
        plot_epsname2 = fullfile(taskDir, sprintf('fixed_ddf_atlas_2x%s_%s_%s.eps', suffix, task, interval_name));
    else
        plot_fname2 = fullfile(taskDir, sprintf('ddf_atlas_2x%s_%s_%s.png', suffix, task, interval_name));
        plot_epsname2 = fullfile(taskDir, sprintf('ddf_atlas_2x%s_%s_%s.eps', suffix, task, interval_name));
    end
    
    saveas(gcf, plot_fname2);
    set(gcf,'renderer','Painters');
    print('-depsc','-tiff','-r300', '-painters', plot_epsname2);
end

