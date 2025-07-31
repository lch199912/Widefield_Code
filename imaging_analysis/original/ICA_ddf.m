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
icaDir = fullfile(fdir, 'ICA_results');
load(fullfile(fdir, "atlas.mat"));
brainmask = ~isnan(brainmask); % convert into boolean
load(fullfile(fdir, 'locaNMF','properties.mat'));
dfDir = fullfile(fdir, 'df_data');
% load(fullfile(dfDir, 'df_data.mat'));
load(fullfile(dfDir, 'FC_data_all.mat'), 'zscore_stacked')

%% 

zscore_1x_all = zscore_stacked(list_1x);
zscore_2x_all = zscore_stacked(list_2x);
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
% if exist('atlas', 'var')
%     atlas = double(atlas);
%     atlas_edge = edge(atlas, 'Sobel');  % extract outline by Sobel
%     
%     % find atlas edge
%     dx = diff(atlas,1,2) ~= 0; 
%     dy = diff(atlas,1,1) ~= 0; 
%     atlas_outline = padarray(dx, [0 1], 0, 'post') | padarray(dy, [1 0], 0, 'post');
% %     atlas_outline = imdilate(atlas_outline, strel('disk', 2));
% end
% imagesc(atlas_outline)

%% 
node_labels = 1:28;
plot_activity_heatmap(zscore_1x, 20, node_labels, 90, 290, [-3, 4], 'jet', 'zscored Heatmap for 1x trial');
set(gcf, 'WindowState','maximized')
saveas(gcf, fullfile(dfDir, 'zscore_heatmap_1x.png')); 
close(gcf);

plot_activity_heatmap(zscore_2x, 20, node_labels, 90, 290, [-3, 4], 'jet', 'zscored Heatmap for 2x trial');
set(gcf, 'WindowState','maximized')
saveas(gcf, fullfile(dfDir, 'zscore_heatmap_2x.png')); 
close(gcf);
disp('df heatmap saved.')

%% 
% save Mean Δ df 
save(fullfile(dfDir, 'FC_ddf.mat'), 'zs_rest_1x', 'zs_rest_2x', 'mean_zs_1x', 'mean_zs_2x', 'mean_dzs_1x', 'mean_dzs_2x', 'calculate_type','resting_time','walking_time');
disp('Mean Delta df calculated and saved.');

