clear;clc;close all;
%% Set parameters
resting_time = [3, 11]; % resting time index with corr frame (-4s ~ -1s)
walking_time = [19, 27; 39, 47; 59, 67]; % time to compare (0~3/5~8/10~13)
% resting_time = convert_sec_to_window([-3.5, -1.5]);
% walking_time = convert_sec_to_window([0, 2; 5.5, 7.5; 10, 12]);
% plot_option = true; % plot data
zscore_option = true; 

idx_name = {'initiation', 'execution', 'completion'};
region_boundaries = [1, 6; 7, 10; 11, 16; 17, 18; 19, 24; 25, 28]; % Each row [start_node, end_node] defines a region
region_boundaries_LR = [1, 3; 4, 5; 6, 8; 9, 9; 10, 12; 13, 14]; % Each row [start_node, end_node] defines a region
node_labels = {'M2', 'M1', 'S1', 'Aud', 'Vis', 'RSC'}; % Labels for nodes (modify as needed)
%% Load PLSR data
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
load(fullfile(fdir,'locaNMF','properties.mat'));
PLSRDir = fullfile(fdir, 'df_data', 'PLSR_all');
if zscore_option == true
    load(fullfile(PLSRDir, 'PLSR_FC_zscored.mat'), 'filtered_data', 'filt_avg1x', 'filt_avg2x')
    subfix = '_zscored';
elseif zscore_option == false
    load(fullfile(PLSRDir, 'PLSR_FC.mat'), 'filtered_data', 'filt_avg1x', 'filt_avg2x')
    subfix = '';
end
filtered_data_1x = filtered_data(list_1x); filtered_data_2x = filtered_data(list_2x);
colormap_list = load('D:\reference\colormap_data.mat');
%% fisher z-transform
filtered_data_1x_z = cell(size(filtered_data_1x)); 
filtered_data_2x_z = cell(size(filtered_data_2x)); 
delta = 1e-5;  % clipping margin

for i = 1:length(filtered_data_1x)
    fc_mat = filtered_data_1x{i};  % [node x node x window]
    fc_clipped = min(max(fc_mat, -1 + delta), 1 - delta);
    fc_z = atanh(fc_clipped); % Fisher z-transform
    filtered_data_1x_z{i} = fc_z; 
end
for i = 1:length(filtered_data_2x)
    fc_mat = filtered_data_2x{i};  % [node x node x window]
    fc_clipped = min(max(fc_mat, -1 + delta), 1 - delta);
    fc_z = atanh(fc_clipped); % Fisher z-transform
    filtered_data_2x_z{i} = fc_z; 
end
%% extract event data

total_idx = 1 + size(walking_time, 1); % resting + walking phases
num_nodes = size(filtered_data_1x{1},1);
num_trial_1x = length(filtered_data_1x); num_trial_2x = length(filtered_data_2x);

fc_event_avg_1x = cell(total_idx, num_trial_1x);
fc_event_avg_2x = cell(total_idx, num_trial_2x);

%1x
for tidx = 1:num_trial_1x
    data = filtered_data_1x{tidx};
    fc_event_avg_1x{1, tidx} = mean(data(:,:,resting_time(1):resting_time(2)), 3);
    for i = 1:size(walking_time,1)
        fc_event_avg_1x{1+i, tidx} = mean(data(:,:,walking_time(i,1):walking_time(i,2)), 3);
    end
end

%2x
for tidx = 1:num_trial_2x
    data = filtered_data_2x{tidx};
    fc_event_avg_2x{1, tidx} = mean(data(:,:,resting_time(1):resting_time(2)), 3);
    for i = 1:size(walking_time,1)
        fc_event_avg_2x{1+i, tidx} = mean(data(:,:,walking_time(i,1):walking_time(i,2)), 3);
    end
end

mean_R_1x = mean_cell_matrices(filtered_data_1x); mean_R_2x = mean_cell_matrices(filtered_data_2x);
%% delta R (event R - resting R)

delta_fc_1x = cell(size(walking_time,1), num_trial_1x);
delta_fc_2x = cell(size(walking_time,1), num_trial_2x);

for tidx = 1:num_trial_1x
    rest = fc_event_avg_1x{1, tidx};
    for e = 1:size(walking_time,1)
        delta_fc_1x{e, tidx} = fc_event_avg_1x{1+e, tidx} - rest;
    end
end
for tidx = 1:num_trial_2x
    rest = fc_event_avg_2x{1, tidx};
    for e = 1:size(walking_time,1)
        delta_fc_2x{e, tidx} = fc_event_avg_2x{1+e, tidx} - rest;
    end
end

% ΔR left-right difference

left_idx = 1:2:num_nodes;
right_idx = 2:2:num_nodes;

fc_event_avg_1x_LR = cell(size(fc_event_avg_1x));
fc_event_avg_2x_LR = cell(size(fc_event_avg_2x));

for e = 1:size(walking_time, 1)+1  % +1 for resting condition
    for t = 1:num_trial_1x
        L = fc_event_avg_1x{e, t}(left_idx, left_idx);
        R = fc_event_avg_1x{e, t}(right_idx, right_idx);
        fc_event_avg_1x_LR{e, t} = R - L;
    end
    for t = 1:num_trial_2x
        L = fc_event_avg_2x{e, t}(left_idx, left_idx);
        R = fc_event_avg_2x{e, t}(right_idx, right_idx);
        fc_event_avg_2x_LR{e, t} = R - L;
    end
end

delta_fc_1x_LR = cell(size(delta_fc_1x));
delta_fc_2x_LR = cell(size(delta_fc_2x));

for e = 1:size(walking_time,1)
    for t = 1:num_trial_1x
        L = delta_fc_1x{e,t}(left_idx, left_idx);
        R = delta_fc_1x{e,t}(right_idx, right_idx);
        delta_fc_1x_LR{e,t} = R - L;
    end
    for t = 1:num_trial_2x
        L = delta_fc_2x{e,t}(left_idx, left_idx);
        R = delta_fc_2x{e,t}(right_idx, right_idx);
        delta_fc_2x_LR{e,t} = R - L;
    end
end
%% average delta R
fc_avg_1x = zeros(num_nodes, num_nodes, size(walking_time,1)+1);
fc_avg_2x = zeros(num_nodes, num_nodes, size(walking_time,1)+1);
fc_avg_1x_LR = zeros(length(left_idx), length(left_idx), size(walking_time,1)+1);
fc_avg_2x_LR = zeros(length(left_idx), length(left_idx), size(walking_time,1)+1);
for eidx = 1:size(walking_time,1)+1
    fc_avg_1x(:,:,eidx) = mean(cat(3,fc_event_avg_1x{eidx,:}), 3);
    fc_avg_2x(:,:,eidx) = mean(cat(3,fc_event_avg_2x{eidx,:}), 3);
    fc_avg_1x_LR(:,:,eidx) = mean(cat(3,fc_event_avg_1x_LR{eidx,:}), 3);
    fc_avg_2x_LR(:,:,eidx) = mean(cat(3,fc_event_avg_2x_LR{eidx,:}), 3);
end

dfc_avg_1x = zeros(num_nodes, num_nodes, size(walking_time,1));
dfc_avg_2x = zeros(num_nodes, num_nodes, size(walking_time,1));
dfc_avg_1x_LR = zeros(num_nodes/2, num_nodes/2, size(walking_time,1));
dfc_avg_2x_LR = zeros(num_nodes/2, num_nodes/2, size(walking_time,1));

for e = 1:size(walking_time,1)
    dfc_avg_1x(:,:,e) = mean(cat(3, delta_fc_1x{e,:}), 3);
    dfc_avg_2x(:,:,e) = mean(cat(3, delta_fc_2x{e,:}), 3);
    dfc_avg_1x_LR(:,:,e) = mean(cat(3, delta_fc_1x_LR{e,:}), 3);
    dfc_avg_2x_LR(:,:,e) = mean(cat(3, delta_fc_2x_LR{e,:}), 3);
end
%% 
% 사용자 입력
unit_list = 1:6;  % or multiple units
data_to_use = mean_R_1x;  % size = [28 x 28 x 73]
save_dir = fullfile(PLSRDir, 'heatmap_single_unit');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

for idx = 1:length(unit_list)
    u = unit_list(idx);  % source unit

    % 데이터: 모든 target unit과의 correlation 값 (시간 축 따라 변화)
    org_data = squeeze(data_to_use(u, :, :));  % [1 x 28 x 73] → [28 x 73]
    
    % Figure
    fig = figure('Visible', 'on');
    set(gcf, 'WindowState','maximized');
    
    imagesc(org_data);  % y = target unit, x = time
    colormap('jet');
    colorbar;
    caxis([-0.3 0.3]);  % 또는 auto 설정
    xlabel('Time window');
    ylabel('Target Unit Index');
    title(sprintf('Correlation: Unit %d vs All Units Over Time', u));
    
    % 유닛 경계선 추가
    draw_region_boundaries(size(org_data, 1));  % y축 기준
    set(gca, 'YDir', 'normal');  % y축 방향 위에서 아래로

    % 저장
    fname = sprintf('Corr_OverTime_Unit%02d_avg1x.png', u);
    saveas(fig, fullfile(save_dir, fname));

    fname_fig = sprintf('Corr_OverTime_Unit%02d_avg1x.fig', u);
    savefig(fig, fullfile(save_dir, fname_fig));
    close(fig);
end


%% 
parameters = struct();
parameters.resting_time = resting_time;
parameters.walking_time = walking_time;
parameters.idx_name = idx_name;
parameters.region_boundaries = region_boundaries;
parameters.region_boundaries_LR = region_boundaries_LR;
parameters.node_labels = node_labels;
% parameters.plot_option = plot_option;

save(fullfile(PLSRDir, ['PLSR_FC_mean', subfix, '.mat']), 'parameters', 'dfc_avg_1x', 'dfc_avg_2x', 'fc_avg_1x','fc_avg_2x');
disp('Averaged R and dR data calculated and saved.');
save(fullfile(PLSRDir, ['PLSR_FC_LR', subfix, '.mat']), 'fc_avg_1x_LR', 'fc_avg_2x_LR', 'dfc_avg_1x_LR','dfc_avg_2x_LR');
disp('Averaged difference of LR for R and dR data saved.');

save(fullfile(PLSRDir, ['PLSR_FC_raw',subfix,'.mat']), 'delta_fc_1x', 'delta_fc_2x', ...
    'delta_fc_1x_LR', 'delta_fc_2x_LR', 'fc_event_avg_1x', 'fc_event_avg_2x', 'filtered_data_1x','filtered_data_2x','filtered_data_1x_z', 'filtered_data_2x_z');
disp('ΔR calculation done')
