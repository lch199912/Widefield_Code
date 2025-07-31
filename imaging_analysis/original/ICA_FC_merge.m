clear; clc; close all;
%% Base parameters
y_limit_dR = [-0.5, 0.5];
y_limit_R = [0, 1];
alpha = 0.1;
plot_option = false;
idx_name = {'initiation', 'execution', 'completion'};
idx_title = [{'resting'}, idx_name];
region_boundaries = [1, 6; 7, 10; 11, 16; 17, 18; 19, 24; 25, 28]; % Each row [start_node, end_node] defines a region
node_labels = {'M2', 'M1', 'S1', 'Aud', 'Vis', 'RSC'}; % Labels for nodes (modify as needed)
%% Choose task
addpath(genpath('E:\Final Codes'));
load('D:\reference\colormap_data.mat');
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

resultDir = fullfile('D:\data analysis\results', 'avg_dR');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
taskDir = fullfile(resultDir, task);
if ~exist(taskDir, 'dir')
    mkdir(taskDir);
end
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
avg_dR_1x = cell(numFolders, 1); 
avg_dR_2x = cell(numFolders, 1); 
resting_1x = cell(numFolders, 1); 
resting_2x = cell(numFolders, 1); 
avg_R_1x = cell(numFolders, 1); 
avg_R_2x = cell(numFolders, 1); 
diff = cell(numFolders, 1);


for i = 1:numFolders
    topFolder = topFolders{i};
    [~, mouseName] = fileparts(topFolder);
    subFolders = dir(topFolder);
    subFolders = subFolders([subFolders.isdir]);
    subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

    [~, idx_recent] = max([subFolders.datenum]);
    recentFolder = fullfile(topFolder, subFolders(idx_recent).name);

    targetFile = fullfile(recentFolder, 'imaging', 'df_data', 'FC_data_mean.mat');
    targetFile2 = fullfile(recentFolder, 'imaging', 'df_data', 'FC_data_all.mat');

    if exist(targetFile, 'file')
        loadedData = load(targetFile);
        data(i).mouseName = mouseName;
        data(i).data = loadedData;

        if isfield(loadedData, 'dfc_avg_1x')
            avg_dR_1x{i} = loadedData.dfc_avg_1x; 
            avg_dR_2x{i} = loadedData.dfc_avg_2x; 
            avg_R_1x{i} = loadedData.fc_avg_1x; 
            avg_R_2x{i} = loadedData.fc_avg_2x;
            avgZ_dR_1x{i} = loadedData.dfc_avg_z_1x; 
            avgZ_dR_2x{i} = loadedData.dfc_avg_z_2x; 
            avgZ_R_1x{i} = loadedData.fc_avg_z_1x; 
            avgZ_R_2x{i} = loadedData.fc_avg_z_2x;
        else
            warning('data not found in file %s', targetFile);
        end
    else
        warning('File does not exist: %s', targetFile);
    end
end

save(fullfile(taskDir, 'raw_data.mat'), 'avg_R_1x', 'avg_R_2x','avg_dR_1x','avg_dR_2x','avgZ_R_1x','avgZ_R_2x','avgZ_dR_1x','avgZ_dR_2x');
%% 
num_animals = numel(avgZ_dR_1x);
num_nodes = size(avgZ_dR_1x{1}, 1);
num_intervals = size(avgZ_dR_1x{1}, 3);

p_values_dR_1x = NaN(num_nodes, num_nodes, num_intervals);
h_values_dR_1x = NaN(num_nodes, num_nodes, num_intervals);
all_p_values_1x = [];
all_indices_1x = [];

% ΔZ 값에 대한 t-test 수행
for eidx = 1:num_intervals
    for i = 1:num_nodes
        for j = 1:num_nodes
            deltaz_all = zeros(num_animals, 1);
            for aidx = 1:num_animals
                deltaz_all(aidx) = avgZ_dR_1x{aidx}(i, j, eidx);  % ΔZ = interval - rest
            end

            [~, p] = ttest(deltaz_all, 0, 'Alpha', alpha);
            p_values_dR_1x(i, j, eidx) = p;

            all_p_values_1x = [all_p_values_1x; p];
            all_indices_1x = [all_indices_1x; i, j, eidx];
        end
    end
end

% FDR 보정
[h_fdr_dR_1x, crit_p_dR_1x] = fdr_bh(all_p_values_1x, alpha);

% 유의한 인덱스 저장
significant_idx_1x = find(h_fdr_dR_1x == 1);
sorted_idx_dR_1x = all_indices_1x(significant_idx_1x, :);
sorted_p_values_dR_1x = all_p_values_1x(significant_idx_1x);

% 유의성 표시 행렬 생성
for k = 1:length(significant_idx_1x)
    i = sorted_idx_dR_1x(k, 1);
    j = sorted_idx_dR_1x(k, 2);
    eidx = sorted_idx_dR_1x(k, 3);
    h_values_dR_1x(i, j, eidx) = 1;
end
%% 
p_values_dR_2x = NaN(num_nodes, num_nodes, num_intervals);
h_values_dR_2x = NaN(num_nodes, num_nodes, num_intervals);
all_p_values_2x = [];
all_indices_2x = [];

% ΔZ 값에 대한 t-test 수행
for eidx = 1:num_intervals
    for i = 1:num_nodes
        for j = 1:num_nodes
            deltaz_all = zeros(num_animals, 1);
            for aidx = 1:num_animals
                deltaz_all(aidx) = avgZ_dR_2x{aidx}(i, j, eidx);
            end

            [~, p] = ttest(deltaz_all, 0, 'Alpha', alpha);
            p_values_dR_2x(i, j, eidx) = p;

            all_p_values_2x = [all_p_values_2x; p];
            all_indices_2x = [all_indices_2x; i, j, eidx];
        end
    end
end

% FDR 보정
[h_fdr_dR_2x, crit_p_dR_2x] = fdr_bh(all_p_values_2x, alpha);

% 유의한 인덱스 저장
significant_idx_2x = find(h_fdr_dR_2x == 1);
sorted_idx_dR_2x = all_indices_2x(significant_idx_2x, :);
sorted_p_values_dR_2x = all_p_values_2x(significant_idx_2x);

% 유의성 표시 행렬 생성
for k = 1:length(significant_idx_2x)
    i = sorted_idx_dR_2x(k, 1);
    j = sorted_idx_dR_2x(k, 2);
    eidx = sorted_idx_dR_2x(k, 3);
    h_values_dR_2x(i, j, eidx) = 1;
end

%% 

save(fullfile(taskDir, sprintf('ttest_dR_%s.mat', task)), ...
    'p_values_dR_1x', 'h_values_dR_1x', 'sorted_idx_dR_1x', 'sorted_p_values_dR_1x', 'crit_p_dR_1x', ...
    'p_values_dR_2x', 'h_values_dR_2x', 'sorted_idx_dR_2x', 'sorted_p_values_dR_2x', 'crit_p_dR_2x');

%% Average data
mean_R_1x = mean_cell_matrices(avg_R_1x);
mean_R_2x = mean_cell_matrices(avg_R_2x);
mean_dR_1x = mean_cell_matrices(avg_dR_1x);
mean_dR_2x = mean_cell_matrices(avg_dR_2x);
save(fullfile(taskDir, sprintf('mean_corr_all_%s.mat', task)), 'mean_R_1x', 'mean_R_2x','mean_dR_1x','mean_dR_2x');
%% 
if plot_option == false
    return;
end
%% Plot heatmap 1x R
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(mean_R_1x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_R_1x(:,:,i);

    title_str = sprintf('Mean R (1x) of %s - %s', task, idx_title{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
    colormap(flipud(slanCM('RdBu')));  
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_R); % Apply same color limits
    axis equal
    
    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end

    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'R_avg_1x_heatmap.png'));
plot_epsname = fullfile(taskDir, 'R_avg_1x_heatmap.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('FC heatmaps for 1x condition saved.');
%% Plot heatmap 2x R
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(mean_R_2x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_R_2x(:,:,i);
    title_str = sprintf('Mean R (2x) of %s - %s', task, idx_title{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
    colormap(flipud(slanCM('RdBu'))); 
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_R); % Apply same color limits
    axis equal

    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end

    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'R_avg_2x_heatmap.png'));
plot_epsname = fullfile(taskDir, 'R_avg_2x_heatmap.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('FC heatmaps for 2x condition saved.');
%% Plot heatmap 1x dR
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(mean_dR_1x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_dR_1x(:,:,i);
    title_str = sprintf('Mean dR (1x) of %s - %s', task, idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
%     colormap(flipud(slanCM('RdBu'))); 
    colormap(custom_RdBu);
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_dR); % Apply same color limits
    axis equal

    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end
    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'dR_avg_1x_heatmap.png'));
plot_epsname = fullfile(taskDir, 'dR_avg_1x_heatmap.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('delta FC heatmaps for 1x condition saved.');

%% Plot heatmap 2x
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(mean_dR_2x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_dR_2x(:,:,i);
    title_str = sprintf('Mean dR (2x) of %s - %s', task, idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
%     colormap(flipud(slanCM('RdBu'))); 
    colormap(custom_RdBu);
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_dR); % Apply same color limits
    axis equal

    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end
    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'dR_avg_2x_heatmap.png'));
plot_epsname = fullfile(taskDir, 'dR_avg_2x_heatmap.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('delta FC heatmaps for 2x condition saved.');
%% Plot heatmap 1x mask
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(mean_dR_1x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_dR_1x(:,:,i);
    mask = h_values_dR_1x(:,:,i);
    data(isnan(mask)) = 0;  % not significant → mask to 0
    title_str = sprintf('Mean dR (1x) of %s - %s', task, idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
%     colormap(flipud(slanCM('RdBu')));
    colormap(custom_RdBu);
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_dR); % Apply same color limits
    axis equal

    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end
    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'dR_avg_1x_heatmap_mask.png'));
plot_epsname = fullfile(taskDir, 'dR_avg_1x_heatmap_mask.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('delta FC heatmaps for 1x condition saved.');

%% Plot heatmap 2x mask
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(mean_dR_2x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_dR_2x(:,:,i);
    mask = h_values_dR_2x(:,:,i);
    data(isnan(mask)) = 0;  % not significant → mask to 0
    title_str = sprintf('Mean dR (2x) of %s - %s', task, idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
%     colormap(flipud(slanCM('RdBu'))); 
    colormap(custom_RdBu);
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_dR); % Apply same color limits
    axis equal

    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end
    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'dR_avg_2x_heatmap_mask.png'));
plot_epsname = fullfile(taskDir, 'dR_avg_2x_heatmap_mask.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('delta FC heatmaps for 2x condition saved.');
%% Plot heatmap 2x-1x
figure;
set(gcf, 'WindowState', 'maximized');
mean_diff = mean_dR_2x - mean_dR_1x;
y_limit_ddR = [-0.3, 0.3];

for i = 1:size(mean_diff , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = mean_diff(:,:,i);
    title_str = sprintf('Mean fraction dR(2x-1x) of %s - %s', task, idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
%     BWP_colormap = customcolormap(linspace(0,1,11), ...
%             {'#523107','#523107','#bf812c','#e2c17e','#f3e9c4','#f6f4f4','#cae9e3','#81cdc1','#379692','#01665e','#003d2e'});
%     colormap(BWP_colormap); 
%     colormap(flipud(slanCM('RdBu')));
    colormap(custom_RdBu);
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(y_limit_ddR); % Apply same color limits
    axis equal

    % Add region boundaries (black lines **between** regions)
    hold on;
    [nrows, ncols] = size(data);
    for r = 1:nrows
        for c = 1:ncols
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.1);  % light gray thin lines
        end
    end
    for r = 1:size(region_boundaries, 1)-1
        boundary = region_boundaries(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end
    % Add node labels
    set(gca, 'XTick', mean(region_boundaries, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(taskDir, 'dR_frac_heatmap.png'));
plot_epsname = fullfile(taskDir, 'dR_frac_heatmap.eps');
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', plot_epsname);
disp('delta FC heatmaps for fraction (2x-1x) condition saved.');