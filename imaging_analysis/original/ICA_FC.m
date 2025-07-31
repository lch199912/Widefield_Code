clear;clc;close all;
%% parameters
window_size = 20; % 1s window (20 frames)
step_size = 5; % 0.25s step (5 frames)

% resting_time = [5, 15]; %  (-3s ~ 0s)
% walking_time = [19, 35; 39, 55; 59, 71]; % (0~5/5~10/10~14)
resting_time = [3, 11]; % resting time index with corr frame (-4s ~ -1s)
walking_time = [19, 27; 39, 47; 59, 67]; % time to compare (0~3/5~8/10~13)
% walking_time = [19, 27; 33, 41; 59, 67]; % time to compare (0~3/3.5~6.5/10~13)
% walking_time = convert_sec_to_window([3.5, 6.5; 5, 7; 5, 8]);
% walking_time = [19, 27; 38, 45; 59, 67];
% resting_time = convert_sec_to_window([-3.5, -1.5]);
% walking_time = convert_sec_to_window([0, 2; 5.5, 7.5; 10, 12]);
plot_option = false; % plot data
raw_option = true; % save raw data

caxis_limits_R = [-1 1]; % Adjust color range for R value plot
caxis_limits_dR = [-0.6 0.6]; % Adjust color range for dR value plot
caxis_limits_R_LR = [-0.3 0.3]; % Adjust color range for R value plot
caxis_limits_dR_LR = [-0.3 0.3]; % Adjust color range for dR value plot

idx_name = {'initiation', 'execution', 'completion'};
region_boundaries = [1, 6; 7, 10; 11, 16; 17, 18; 19, 24; 25, 28]; % Each row [start_node, end_node] defines a region
region_boundaries_LR = [1, 3; 4, 5; 6, 8; 9, 9; 10, 12; 13, 14]; % Each row [start_node, end_node] defines a region
node_labels = {'M2', 'M1', 'S1', 'Aud', 'Vis', 'RSC'}; % Labels for nodes (modify as needed)
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
icaDir = fullfile(fdir, 'ICA_results');
load(fullfile(fdir, "atlas.mat"));
brainmask = ~isnan(brainmask); % convert into boolean
load(fullfile(fdir,'locaNMF','properties.mat'));
dfDir = fullfile(fdir, 'df_data');
load(fullfile(dfDir, 'df_data.mat'));
load('D:\reference\BWR_colormap.mat');
load(fullfile(dfDir, 'raw_data.mat'),'node_trial');
%% Calculate z-score data for each trial
num_nodes = size(node_trial, 1);
num_trial = size(node_trial, 2);
num_frames = size(node_trial{1,1}, 1);
zscore_all = cell(size(node_trial));

for tidx = 1:num_trial
    node_data = cell2mat(node_trial(:, tidx)');
    zscore_data = zscore(node_data, 0, 1);
    for nidx = 1:num_nodes
        zscore_all{nidx, tidx} = zscore_data(:, nidx);
    end
end
% Preallocate new z-score cell array
zscore_stacked = cell(1, num_trial);

for tidx = 1:num_trial
    % Stack all nodes for the current trial into a matrix of size frame x node
    zscore_stacked{1, tidx} = cell2mat(reshape(zscore_all(:, tidx), 1, []));
end

%% Functional Connectivity
num_windows = floor((num_frames - window_size) / step_size) + 1;
fc_all = cell(size(zscore_stacked));
fc_all_z = cell(size(zscore_stacked));

for tidx = 1:num_trial
    curr_trial = zscore_stacked{tidx};
    fc = zeros(num_nodes, num_nodes, num_windows);
    for w = 1:num_windows
        start_idx = (w - 1) * step_size + 1;
        end_idx = start_idx + window_size - 1;
        
        window_data = curr_trial(start_idx:end_idx, :);
        fc(:,:,w) = corr(window_data, 'Type', 'Pearson', 'Rows', 'pairwise');
    end
    fc_all{tidx} = fc;

    % fisher z-transform
    delta = 1e-5;
    R_clipped = min(max(fc, -1 + delta), 1 - delta);
    Z = atanh(R_clipped);

    fc_all_z{tidx} = Z;
end

%% extract event data for fc_all
total_idx = 1 + size(walking_time, 1);
fc_event_avg = cell(total_idx, num_trial);
for tidx = 1:num_trial
    curr_trial = fc_all{tidx};
    resting = curr_trial(:,:,resting_time(1):resting_time(2));
    resting_mean = mean(resting, 3, 'omitnan');
    fc_event_avg{1, tidx} = resting_mean;
    
    for event_idx = 1:size(walking_time, 1)
        event = curr_trial(:,:,walking_time(event_idx,1):walking_time(event_idx,2));
        event_mean = mean(event, 3, 'omitnan');
        fc_event_avg{event_idx+1, tidx} = event_mean;
    end
end

%% delta (event - resting) data for fc_all
delta_fc = cell(size(walking_time,1), num_trial);
for tidx = 1:num_trial
    curr_resting = fc_event_avg{1, tidx};
    for event_idx = 1:size(walking_time,1)
        event = fc_event_avg{1+event_idx, tidx};
        devent = event - curr_resting;
        delta_fc{event_idx, tidx} = devent;
    end
end

%% extract event data for fc_all_z
fc_event_avg_z = cell(total_idx, num_trial);
for tidx = 1:num_trial
    curr_trial = fc_all_z{tidx};
    resting = curr_trial(:,:,resting_time(1):resting_time(2));
    resting_mean = mean(resting, 3, 'omitnan');
    fc_event_avg_z{1, tidx} = resting_mean;
    
    for event_idx = 1:size(walking_time, 1)
        event = curr_trial(:,:,walking_time(event_idx,1):walking_time(event_idx,2));
        event_mean = mean(event, 3, 'omitnan');
        fc_event_avg_z{event_idx+1, tidx} = event_mean;
    end
end

%% delta (event - resting) data for fc_all_z
delta_fc_z = cell(size(walking_time,1), num_trial);
for tidx = 1:num_trial
    curr_resting = fc_event_avg_z{1, tidx};
    for event_idx = 1:size(walking_time,1)
        event = fc_event_avg_z{1+event_idx, tidx};
        devent = event - curr_resting;
        delta_fc_z{event_idx, tidx} = devent;
    end
end


%% distinguish 1x and 2x trial, then average

fc_event_avg_1x = fc_event_avg(:,list_1x);
fc_event_avg_2x = fc_event_avg(:,list_2x);
delta_fc_1x = delta_fc(:,list_1x);
delta_fc_2x = delta_fc(:,list_2x);
fc_avg_1x = zeros(num_nodes, num_nodes, size(walking_time,1)+1);
fc_avg_2x = zeros(num_nodes, num_nodes, size(walking_time,1)+1);
dfc_avg_1x = zeros(num_nodes, num_nodes, size(walking_time,1));
dfc_avg_2x = zeros(num_nodes, num_nodes, size(walking_time,1));

for eidx = 1:size(walking_time,1)+1
    fc_avg_1x(:,:,eidx) = mean(cat(3,fc_event_avg_1x{eidx,:}), 3);
    fc_avg_2x(:,:,eidx) = mean(cat(3,fc_event_avg_2x{eidx,:}), 3);
end
for eidx = 1:size(walking_time,1)
    dfc_avg_1x(:,:,eidx) = mean(cat(3,delta_fc_1x{eidx,1}), 3);
    dfc_avg_2x(:,:,eidx) = mean(cat(3,delta_fc_2x{eidx,1}), 3);
end
%% distinguish 1x and 2x trial, then average (Z-scored data)

fc_event_avg_z_1x = fc_event_avg_z(:,list_1x);
fc_event_avg_z_2x = fc_event_avg_z(:,list_2x);
delta_fc_z_1x = delta_fc_z(:,list_1x);
delta_fc_z_2x = delta_fc_z(:,list_2x);

fc_avg_z_1x = zeros(num_nodes, num_nodes, size(walking_time,1)+1);
fc_avg_z_2x = zeros(num_nodes, num_nodes, size(walking_time,1)+1);
dfc_avg_z_1x = zeros(num_nodes, num_nodes, size(walking_time,1));
dfc_avg_z_2x = zeros(num_nodes, num_nodes, size(walking_time,1));

for eidx = 1:size(walking_time,1)+1
    fc_avg_z_1x(:,:,eidx) = mean(cat(3,fc_event_avg_z_1x{eidx,:}), 3);
    fc_avg_z_2x(:,:,eidx) = mean(cat(3,fc_event_avg_z_2x{eidx,:}), 3);
end
for eidx = 1:size(walking_time,1)
    dfc_avg_z_1x(:,:,eidx) = mean(cat(3,delta_fc_z_1x{eidx,:}), 3);
    dfc_avg_z_2x(:,:,eidx) = mean(cat(3,delta_fc_z_2x{eidx,:}), 3);
end

%% LR comparison
num_nodes_single = num_nodes/2;
left_hemi_idx = 1:2:num_nodes; right_hemi_idx = 2:2:num_nodes;

fc_event_avg_1x_L = cell(size(fc_event_avg_1x)); fc_event_avg_1x_R = cell(size(fc_event_avg_1x));
fc_event_avg_2x_L = cell(size(fc_event_avg_2x)); fc_event_avg_2x_R = cell(size(fc_event_avg_2x));
for eidx = 1:size(walking_time,1)+1
    for tidx = 1:length(list_1x)
        curr_1x = fc_event_avg_1x{eidx, tidx};
        curr_1x_L = curr_1x(left_hemi_idx, left_hemi_idx);
        curr_1x_R = curr_1x(right_hemi_idx, right_hemi_idx);
        fc_event_avg_1x_L{eidx, tidx} = curr_1x_L; fc_event_avg_1x_R{eidx, tidx} = curr_1x_R;
    end
    for tidx = 1:length(list_2x)
        curr_2x = fc_event_avg_2x{eidx, tidx};
        curr_2x_L = curr_2x(left_hemi_idx, left_hemi_idx);
        curr_2x_R = curr_2x(right_hemi_idx, right_hemi_idx);
        fc_event_avg_2x_L{eidx, tidx} = curr_2x_L; fc_event_avg_2x_R{eidx, tidx} = curr_2x_R;
    end
end

%% LR comparison for delta_fc
delta_fc_1x_L = cell(size(delta_fc_1x)); delta_fc_1x_R = cell(size(delta_fc_1x));
delta_fc_2x_L = cell(size(delta_fc_2x)); delta_fc_2x_R = cell(size(delta_fc_2x));

for eidx = 1:size(walking_time,1)
    for tidx = 1:length(list_1x)
        curr_1x = delta_fc_1x{eidx, tidx};
        curr_1x_L = curr_1x(left_hemi_idx, left_hemi_idx);
        curr_1x_R = curr_1x(right_hemi_idx, right_hemi_idx);
        delta_fc_1x_L{eidx, tidx} = curr_1x_L; delta_fc_1x_R{eidx, tidx} = curr_1x_R;
    end
    for tidx = 1:length(list_2x)
        curr_2x = delta_fc_2x{eidx, tidx};
        curr_2x_L = curr_2x(left_hemi_idx, left_hemi_idx);
        curr_2x_R = curr_2x(right_hemi_idx, right_hemi_idx);
        delta_fc_2x_L{eidx, tidx} = curr_2x_L; delta_fc_2x_R{eidx, tidx} = curr_2x_R;
    end
end
%% LR difference calculation
% Compute LR differences for fc_event_avg
fc_event_avg_1x_LR = cell(size(fc_event_avg_1x_L));
fc_event_avg_2x_LR = cell(size(fc_event_avg_2x_L));

for eidx = 1:size(walking_time,1)+1
    for tidx = 1:length(list_1x)
        fc_event_avg_1x_LR{eidx, tidx} = fc_event_avg_1x_R{eidx, tidx} - fc_event_avg_1x_L{eidx, tidx};
    end
    for tidx = 1:length(list_2x)
        fc_event_avg_2x_LR{eidx, tidx} = fc_event_avg_2x_R{eidx, tidx} - fc_event_avg_2x_L{eidx, tidx};
    end
end

% Compute LR differences for delta_fc
delta_fc_1x_LR = cell(size(delta_fc_1x_L));
delta_fc_2x_LR = cell(size(delta_fc_2x_L));

for eidx = 1:size(walking_time,1)
    for tidx = 1:length(list_1x)
        delta_fc_1x_LR{eidx, tidx} = delta_fc_1x_R{eidx, tidx} - delta_fc_1x_L{eidx, tidx};
    end
    for tidx = 1:length(list_2x)
        delta_fc_2x_LR{eidx, tidx} = delta_fc_2x_R{eidx, tidx} - delta_fc_2x_L{eidx, tidx};
    end
end
%% Average LR difference data

% Initialize mean matrices for LR difference
fc_avg_1x_LR = zeros(num_nodes_single, num_nodes_single, size(walking_time,1)+1);
fc_avg_2x_LR = zeros(num_nodes_single, num_nodes_single, size(walking_time,1)+1);
dfc_avg_1x_LR = zeros(num_nodes_single, num_nodes_single, size(walking_time,1));
dfc_avg_2x_LR = zeros(num_nodes_single, num_nodes_single, size(walking_time,1));

% Compute mean for fc_event_avg_LR (R - L)
for eidx = 1:size(walking_time,1)+1
    fc_avg_1x_LR(:,:,eidx) = mean(cat(3,fc_event_avg_1x_LR{eidx,:}), 3);
    fc_avg_2x_LR(:,:,eidx) = mean(cat(3,fc_event_avg_2x_LR{eidx,:}), 3);
end

% Compute mean for delta_fc_LR (R - L)
for eidx = 1:size(walking_time,1)
    dfc_avg_1x_LR(:,:,eidx) = mean(cat(3,delta_fc_1x_LR{eidx,:}), 3);
    dfc_avg_2x_LR(:,:,eidx) = mean(cat(3,delta_fc_2x_LR{eidx,:}), 3);
end

%% save data
parameters = struct();
parameters.window_size = window_size;
parameters.step_size = step_size;
parameters.resting_time = resting_time;
parameters.walking_time = walking_time;
parameters.idx_name = idx_name;
parameters.region_boundaries = region_boundaries;
parameters.region_boundaries_LR = region_boundaries_LR;
parameters.node_labels = node_labels;
parameters.plot_option = plot_option;
parameters.raw_option = raw_option;

if raw_option == true
    save(fullfile(dfDir, 'FC_data_all.mat'), 'fc_all_z', 'fc_all', 'zscore_stacked', 'fc_event_avg', 'delta_fc');
    disp('Save zscore data and fc data for all trial...');
    save(fullfile(dfDir, 'FC_data_all_LR.mat'), 'fc_event_avg_1x_L', 'fc_event_avg_1x_R', 'fc_event_avg_2x_L', 'fc_event_avg_2x_R', ...
        'delta_fc_1x_L','delta_fc_1x_R','delta_fc_2x_L', 'delta_fc_2x_R');
    disp('Save whole LR difference data...');
end

save(fullfile(dfDir, 'FC_data_mean.mat'), 'parameters', 'dfc_avg_1x', 'dfc_avg_2x', 'fc_avg_1x','fc_avg_2x', 'fc_avg_z_1x', 'fc_avg_z_2x','dfc_avg_z_1x','dfc_avg_z_2x');
disp('Averaged R and dR data calculated and saved.');
save(fullfile(dfDir, 'FC_data_mean_LR.mat'), 'fc_avg_1x_LR', 'fc_avg_2x_LR', 'dfc_avg_1x_LR','dfc_avg_2x_LR');
disp('Averaged difference of LR for R and dR data saved.');
%% 
if plot_option == false
    return;
end

%% Plot data
figure;
set(gcf, 'WindowState', 'maximized');

idx_title = [{'resting'}, idx_name];
for i = 1:size(fc_avg_1x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = fc_avg_1x(:,:,i);
    title_str = sprintf('Mean R (1x) - %s', idx_title{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
    colormap(BWR_colormap); 
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(caxis_limits_R); % Apply same color limits

    % Add region boundaries (black lines **between** regions)
    hold on;
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
saveas(gcf, fullfile(dfDir, 'fc_avg_1x_heatmap.png'));

disp('FC heatmaps for 1x condition saved.');

%% 
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(dfc_avg_1x , 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = dfc_avg_1x(:,:,i);
    title_str = sprintf('Mean dR (1x) - %s', idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
    colormap(BWR_colormap); 
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(caxis_limits_dR); % Apply same color limits

    % Add region boundaries (black lines **between** regions)
    hold on;
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
saveas(gcf, fullfile(dfDir, 'dfc_avg_1x_heatmap.png'));

disp('delta FC heatmaps for 1x condition saved.');

%% Plot LR difference (R - L) for fc_avg_1x_LR
figure;
set(gcf, 'WindowState', 'maximized');

idx_title = [{'resting'}, idx_name];
for i = 1:size(fc_avg_1x_LR, 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = fc_avg_1x_LR(:,:,i);
    title_str = sprintf('Mean R (1x) - LR Diff - %s', idx_title{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
    colormap(BWR_colormap); 
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(caxis_limits_R_LR); % Apply same color limits as other FC plots

    % Add region boundaries (black lines **between** regions)
    hold on;
    for r = 1:size(region_boundaries_LR, 1)-1
        boundary = region_boundaries_LR(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end

    % Add node labels
    set(gca, 'XTick', mean(region_boundaries_LR, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries_LR, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(dfDir, 'fc_avg_1x_LR_heatmap.png'));

disp('FC heatmaps for 1x condition (LR difference) saved.');

%% Plot LR difference (R - L) for delta_fc_1x_LR
figure;
set(gcf, 'WindowState', 'maximized');

for i = 1:size(dfc_avg_1x_LR, 3)
    subplot(2,2,i);
    
    % Select data for current subplot
    data = dfc_avg_1x_LR(:,:,i);
    title_str = sprintf('Mean dR (1x) - LR Diff - %s', idx_name{i});
    
    % Set diagonal (self-correlation) to BLACK (-1 ensures it is mapped to black in colormap)
    data(logical(eye(size(data)))) = -1;
    
    % Display heatmap
    imagesc(data);
    colormap(BWR_colormap); 
    colorbar;
    title(title_str);
    xlabel('Node Index');
    ylabel('Node Index');
    caxis(caxis_limits_dR_LR); % Apply same color limits as other delta R plots

    % Add region boundaries (black lines **between** regions)
    hold on;
    for r = 1:size(region_boundaries_LR, 1)-1
        boundary = region_boundaries_LR(r,2) + 0.5;
        plot([0.5, size(data,1)+0.5], [boundary, boundary], 'k', 'LineWidth', 2); % Horizontal
        plot([boundary, boundary], [0.5, size(data,1)+0.5], 'k', 'LineWidth', 2); % Vertical
    end

    % Add node labels
    set(gca, 'XTick', mean(region_boundaries_LR, 2), 'XTickLabel', node_labels, 'FontSize', 10);
    set(gca, 'YTick', mean(region_boundaries_LR, 2), 'YTickLabel', node_labels, 'FontSize', 10);
    
    hold off;
end

% Save Figure
saveas(gcf, fullfile(dfDir, 'dfc_avg_1x_LR_heatmap.png'));

disp('delta FC heatmaps for 1x condition (LR difference) saved.');
