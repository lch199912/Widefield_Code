clear;clc;close all;
%% Choose task
save_eps = true; option_plot = false; zscore_option = true; pupil_scaling_option = 'zscore'; % none / norm / zscore / both
addpath(genpath('E:\Final Codes'));
data_list = {'z-score', 'FC'};
data_choice = questdlg('Select Task:', 'Task Selection', ...
    'z-score', 'FC', 'z-score'); % default == 'treadmill'
if isempty(data_choice)
    disp('Data type selection canceled. Exiting...');
    return;
end
data_map = containers.Map({'z-score', 'FC'}, {'z-score', 'FC'});
data_type = data_map(data_choice);

disp(['Selected Task: ', data_type]);
%% Load behavior and brain data (z-score activity or FC data)
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
load(fullfile(fdir, 'locaNMF','properties.mat'));
dataDir = fullfile(fdir, 'df_data');
if strcmpi(data_type, 'z-score')
    load(fullfile(dataDir, 'FC_data_all.mat'), 'zscore_stacked');
    data = zscore_stacked;
    clearvars zscore_stacked
elseif strcmpi(data_type, 'FC')
    load(fullfile(dataDir, 'FC_data_all.mat'), 'fc_all');
    data = fc_all;
    clearvars fc_all 
    n_trials = length(data);
    n_nodes = size(data{1}, 1);
    idx_lower = find(tril(true(n_nodes), -1));
    processed_data = cell(1, n_trials);
    for i = 1:n_trials
        trial_data = data{i};  % [n_nodes x n_nodes x n_frames]
        n_frames = size(trial_data, 3);
        lower_values = zeros(n_frames, length(idx_lower));
        for f = 1:n_frames
            fc_mat = trial_data(:, :, f);  % 1 frame의 FC
            lower_values(f, :) = fc_mat(idx_lower)'; 
        end
        processed_data{i} = lower_values;  % [n_frames x n_lower_entries]
    end
    data_org = data;
    data_org_1x = data(:,list_1x); data_org_2x = data(:,list_2x);
    data = processed_data;
end
if ~exist('motorSpeed_img', 'var')
    motorSpeed_img = motorSpeed;
end
if length(data) > length(motorSpeed_img)
    data = data(1:length(motorSpeed_img));
end

behavior_dir = fullfile(analysis_dir, mouse, day, 'behavior');
face_file = load(fullfile(behavior_dir, 'face_data.mat'), 'raw_data'); 
face_data = face_file.raw_data.pupil_trials; clearvars face_file
load(fullfile(behavior_dir, 'hindpaw_swing_stand_data.mat'),'smoothed_hindpaw_x_all','smoothed_hindpaw_y_all');
if zscore_option == true
    smoothed_hindpaw_x_all_org = smoothed_hindpaw_x_all;
    smoothed_hindpaw_x_all = zscore(smoothed_hindpaw_x_all);
    smoothed_hindpaw_y_all_org = smoothed_hindpaw_y_all;
    smoothed_hindpaw_y_all = zscore(smoothed_hindpaw_y_all);
end
hindpaw_x = smoothed_hindpaw_x_all; hindpaw_y = smoothed_hindpaw_y_all;
if length(face_data) ~= length(motorSpeed); face_data = face_data(:, 1:length(motorSpeed)); end
if size(hindpaw_x,1) ~= length(motorSpeed); hindpaw_x = hindpaw_x(1:length(motorSpeed),:); hindpaw_y = hindpaw_y(1:length(motorSpeed),:);end
% if length(face_data) >= 100; face_data = face_data(:, 1:100); end
% if size(hindpaw_x,1) >= 100; hindpaw_x = hindpaw_x(1:100,:); hindpaw_y = hindpaw_y(1:100,:);end
if ~isnan(dropTrials)
    face_data(dropTrials) = []; hindpaw_x(dropTrials,:) = []; hindpaw_y(dropTrials,:) = [];
end
%% set epoch range

epoch_range = [];  % frame index (based on brain data). if empty use whole frames.
use_epoch = ~isempty(epoch_range);
mean_option = false;             % true = use epoch mean as input
%% convert brain data
num_trial = numel(data);

if use_epoch
    if mean_option
        brain_data = mean_cell_matrices(cellfun(@(x) x(epoch_range,:), data, 'UniformOutput', false)');
    else
        brain_data = cell2mat(cellfun(@(x) x(epoch_range,:), data, 'UniformOutput', false)');
    end
else
    if mean_option
        brain_data = cell2mat(cellfun(@(x) mean(x, 1), data, 'UniformOutput', false)');
    else
        brain_data = cell2mat(data');
    end
end
%% convert pupil data

% 1. convert all data (n_trial x 800)
pupil_all = cell2mat(face_data');
if any(isnan(pupil_all(:)))
    X = pupil_all;
    for col = 1:size(X,2)
        this_col = X(:,col);
        this_col = fillmissing(this_col, 'linear');
%         this_col = fillmissing(this_col, 'nearest');
        X(:,col) = this_col;
    end
    pupil_all = X;
end

switch lower(pupil_scaling_option)
    case 'none'
        % No preprocessing
        disp('No scaling applied to pupil data.');

    case 'norm'
        % Normalize to [0, 1] range by global max
        max_val = max(pupil_all);
        pupil_all = pupil_all / max_val;
        disp('Applied global max normalization to pupil data.');

    case 'zscore'
        % Global z-score across all trials
        z_all_pupil = zscore(pupil_all);
        pupil_all = reshape(z_all_pupil(1:numel(pupil_all)), size(pupil_all));
        disp('Applied global z-score to pupil data.');

    case 'both'
        % First normalize, then z-score
        normed_all_pupil = pupil_all / max(pupil_all);
        z_normed_all_pupil = zscore(normed_all_pupil);
        pupil_all = reshape(z_normed_all_pupil(1:numel(pupil_all)), size(pupil_all));
        disp('Applied both normalization and z-score to pupil data.');

    otherwise
        error('Invalid pupil_scaling_option. Choose from: none, norm, zscore, both');
end

% 3. remove 20 frames → [n_trial x 760]
pupil_cut = pupil_all(:, 21:780);

% 4. Downsample (e.g. 2 frame interval → 380 frame)
if strcmpi(data_type, 'z-score')
    pupil_downsampled = pupil_cut(:, 1:2:end);
elseif strcmpi(data_type, 'FC')
    pupil_downsampled = resample_frames(pupil_cut, n_frames);

end
if use_epoch
    pupil_downsampled = pupil_downsampled(:, epoch_range);
end

% 5. reshape to vector for PLSR
if mean_option
    pupil_input = mean(pupil_downsampled, 1);
else
    pupil_input = reshape(pupil_downsampled', 1, []);
end

%% convert hindpaw data

% 1. calculate dx, dy, speed = sqrt(dx^2 + dy^2)
dx = diff(hindpaw_x, 1, 2);
dy = diff(hindpaw_y, 1, 2);
speed = sqrt(dx.^2 + dy.^2);
hindpaw_speed = [zeros(size(speed,1),1), speed];

if zscore_option == true
    all_speed = hindpaw_speed(:);
    z_all_speed = zscore(all_speed);
    hindpaw_speed = reshape(z_all_speed(1:numel(hindpaw_speed)), size(hindpaw_speed));
end

% 2. remove 20 frames (21~780 → 760 frame)
hindpaw_x_cut = hindpaw_x(:, 21:780);
hindpaw_y_cut = hindpaw_y(:, 21:780);
hindpaw_speed_cut = hindpaw_speed(:, 21:780);

% 3. Downsampling
if strcmpi(data_type, 'z-score')
    hindpaw_x_downsampled = hindpaw_x_cut(:, 1:2:end);
    hindpaw_y_downsampled = hindpaw_y_cut(:, 1:2:end);
    hindpaw_speed_downsampled = hindpaw_speed_cut(:, 1:2:end);
elseif strcmpi(data_type, 'FC')
    hindpaw_x_downsampled = resample_frames(hindpaw_x_cut, n_frames);
    hindpaw_y_downsampled = resample_frames(hindpaw_y_cut, n_frames);
    hindpaw_speed_downsampled = resample_frames(hindpaw_speed_cut, n_frames);
end

if use_epoch
    hindpaw_x_downsampled = hindpaw_x_downsampled(:, epoch_range);
    hindpaw_y_downsampled = hindpaw_y_downsampled(:, epoch_range);
    hindpaw_speed_downsampled = hindpaw_speed_downsampled(:, epoch_range);
end

% 4. vectorize
if mean_option
    hindpaw_x_input = mean(hindpaw_x_downsampled, 1);
    hindpaw_y_input = mean(hindpaw_y_downsampled, 1);
    hindpaw_speed_input = mean(hindpaw_speed_downsampled, 1);
else
    hindpaw_x_input = reshape(hindpaw_x_downsampled', 1, []);
    hindpaw_y_input = reshape(hindpaw_y_downsampled', 1, []);
    hindpaw_speed_input = reshape(hindpaw_speed_downsampled', 1, []);
end

%% behavior matrix
behavior_input = [pupil_input', hindpaw_x_input', hindpaw_y_input', hindpaw_speed_input'];    % [trial*frame x 3]
if any(isnan(behavior_input(:)))
    % Copy to preserve original
    X = behavior_input;
    % Iterate over columns
    for col = 1:size(X,2)
        this_col = X(:,col);
        this_col = fillmissing(this_col, 'linear');
%         this_col = fillmissing(this_col, 'nearest');
        X(:,col) = this_col;
    end
    % Update the variable
    behavior_input = X;
end
%% set parameters for PLSR

parameters.dataset.explanatoryVariables = behavior_input;    % size: [N_obs x N_vars]
parameters.dataset.responseVariables = brain_data;           % size: [N_obs x N_response_vars]

parameters.comparison_type = 'continuous';                      % or 'categorical'

parameters.findBestNComponents = true;
parameters.ncomponents_max = 10;
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.stratify = false;
parameters.contiguous_partitions = true;

parameters.permutationGeneration = true;
parameters.n_permutations = 1000;
parameters.useBootstrapping = false;

% optional
parameters.useKneePoint = false;                                % or useDifference, minDifference 등
parameters.run_with_max_components = false;
%% 
results = PLSR_forRunAnalysis(parameters);
%% 
X1 = behavior_input;                      % [9000 × 3]
X1_aug = [ones(size(X1,1),1), X1];             % [9000 × 4]
BETA1 = results.results.BETA;                % [4 × 28]
Y1_predicted = X1_aug * BETA1;  
Y1 = brain_data; 

% % 1. Mean Squared Error (MSE)
% mse = mean((Y - Y_predicted).^2, 1);  % [1 x n_response_vars]
% disp('MSE per output variable:');
% disp(mse);
% 
% % 2. R² (coefficient of determination)
% SS_res = sum((Y - Y_predicted).^2);
% SS_tot = sum((Y - mean(Y)).^2);
% R2 = 1 - (SS_res ./ SS_tot);  % [1 x n_response_vars]
% disp('R² per output variable:');
% disp(R2);
% 
% % 3. Correlation between predicted and actual
% corrs = diag(corr(Y, Y_predicted));  % [n_response_vars x 1]
% disp('Correlation per output variable:');
% disp(corrs);

%% filtered Y
filtered_Y1 = Y1 - Y1_predicted; 
if isempty(epoch_range)
    frames_per_trial = cellfun(@(x) size(x,1), data);
    frames_per_trial = frames_per_trial(1);
else
    frames_per_trial = length(epoch_range);  
end

%% 
if strcmpi(data_type, 'z-score')
    % Split filtered_Y back into cell array matching original data_1x
    filtered_data = cell(1, num_trial);
    start_idx = 1;
    for i = 1:num_trial
        f = frames_per_trial;
        filtered_data{i} = filtered_Y1(start_idx:start_idx+f-1, :);  % [frames x sources]
        start_idx = start_idx + f;
    end
elseif strcmpi(data_type, 'FC')
    % FC lower triangle index
    lower_idx = find(tril(true(n_nodes), -1));
    n_sources = length(lower_idx);  
    
    % 1x 
    filtered_data = cell(1, num_trial);
    start_idx = 1;
    for i = 1:num_trial
        trial_data = filtered_Y1(start_idx:start_idx + frames_per_trial - 1, :);  % [73 x 378]
        temp = zeros(frames_per_trial, n_nodes, n_nodes);
    
        for t = 1:frames_per_trial
            mat = zeros(n_nodes);  % [28 x 28]
            mat(lower_idx) = trial_data(t, :);  % lower triangle
            mat = mat + mat';  % symmetric
            temp(t, :, :) = mat;
        end
    
        % Transpose to [nodes x nodes x frames]
        filtered_data{i} = permute(temp, [2, 3, 1]);
        start_idx = start_idx + frames_per_trial;
    end
end
disp(['Filtered data split into ' num2str(num_trial) ' trials.']);
%% 

resultDir = fullfile(dataDir, 'PLSR_all');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
if strcmpi(data_type, 'z-score')
    data_1x = data(:, list_1x); data_2x = data(:, list_2x);
    avg1x = mean_cell_matrices(data_1x);avg2x = mean_cell_matrices(data_2x);
elseif strcmpi(data_type, 'FC')
    avg1x = mean_cell_matrices(data_org_1x);avg2x = mean_cell_matrices(data_org_2x);
end
filt_avg1x = mean_cell_matrices(filtered_data(:, list_1x)); filt_avg2x = mean_cell_matrices(filtered_data(:, list_2x)); 
%% 
% return;
if zscore_option == true
    subfix = '_zscored';
else
    subfix = '';
end
parameters = rmfield(parameters, 'dataset');
save(fullfile(resultDir, ['PLSR_', data_type, subfix, '.mat']), 'parameters', 'results', 'filtered_data', 'filt_avg1x', 'filt_avg2x');
save(fullfile(resultDir, ['PLSR_raw_', data_type, subfix, '.mat']), 'X1_aug','BETA1', 'Y1');
%% 
% winopen(resultDir)
% return;
%% 
if strcmpi(data_type, 'z-score')
    figure;
    set(gcf, 'WindowState','maximized');
    org_data = avg1x(11:370, :);
    filt_data = filt_avg1x(11:370, :);
    for i = 1:size(filt_data, 2)
        filt_data(:,i) = smooth(filt_data(:,i));
    end

    subplot(2,1,1)
    imagesc(org_data')
    title('original data')
    colorbar;
    colormap('jet');
    caxis([-3, 3])
    xline(80.5, '--k', 'LineWidth', 1.5);
    xline(280.5, '--k', 'LineWidth', 1.5);
    draw_region_boundaries(size(org_data, 2));
    
    subplot(2,1,2)
    imagesc(filt_data')
    title('PLSR applied data')
    colorbar;
    colormap('jet');
    caxis([-3, 3])
    xline(80.5, '--k', 'LineWidth', 1.5);
    xline(280.5, '--k', 'LineWidth', 1.5);
    draw_region_boundaries(size(filt_data, 2));
    
    saveas(gcf, fullfile(resultDir, ['PLSR_applied_z-score_1x', subfix,'.png']));
    if save_eps
        plot_epsname = fullfile(resultDir, ['PLSR_applied_z-score_1x', subfix,'.eps']);
        set(gcf,'renderer','Painters');
        print('-depsc','-tiff','-r300', '-painters', plot_epsname);
    end

    figure;
    set(gcf, 'WindowState','maximized');
    org_data = avg2x(11:370, :);
    filt_data = filt_avg2x(11:370, :);
    for i = 1:size(filt_data, 2)
        filt_data(:,i) = smooth(filt_data(:,i));
    end
    subplot(2,1,1)
    imagesc(org_data')
    title('original data')
    colorbar;
    colormap('jet');
    caxis([-3, 3])
    xline(80.5, '--k', 'LineWidth', 1.5);
    xline(280.5, '--k', 'LineWidth', 1.5);
    draw_region_boundaries(size(org_data, 2));

    subplot(2,1,2)
    imagesc(filt_data')
    title('PLSR applied data')
    colorbar;
    colormap('jet');
    caxis([-3, 3])
    xline(80.5, '--k', 'LineWidth', 1.5);
    xline(280.5, '--k', 'LineWidth', 1.5);
    draw_region_boundaries(size(filt_data, 2));

    saveas(gcf, fullfile(resultDir, ['PLSR_applied_z-score_2x', subfix,'.png']));
    if save_eps
        plot_epsname = fullfile(resultDir, ['PLSR_applied_z-score_2x', subfix,'.eps']);
        set(gcf,'renderer','Painters');
        print('-depsc','-tiff','-r300', '-painters', plot_epsname);
    end
elseif strcmpi(data_type, 'FC')
    
end

%% Optional Plot 
if option_plot == false
    return;
end

%% 

randi_trial = randi(size(list_1x,2)); randi_trial = 1;
plot_range = [41,760];

pupil_data = smooth(pupil_all_1x(randi_trial, plot_range(1):plot_range(2)), 20);
pawX_data = smooth(hindpaw_x(randi_trial, plot_range(1):plot_range(2)));
pawY_data = smooth(hindpaw_y(randi_trial, plot_range(1):plot_range(2)));
pawVel_data = smooth(hindpaw_speed(randi_trial, plot_range(1):plot_range(2)),20);
%% 

figure;
subplot(4,1,1)
% plot(pawX_data - pawX_data(1));
plot(pawX_data);
xlim([0, 800])
ylim([-2.5, 2.5]);
% ylim([20, 50]);
subplot(4,1,2)
% plot(pawY_data - pawY_data(1));
plot(pawY_data);
xlim([0, 800])
ylim([-3, 2]);
% ylim([-1, 3]);
subplot(4,1,3)
plot(pawVel_data);
xlim([0, 800])
ylim([-2, 2]);
% ylim([-0.5, 1.5]);
subplot(4,1,4)
plot(pupil_data);
xlim([0, 800])
ylim([-2, 2]);
% ylim([1.5*10^4, 3*10^4]);
set(gcf,'renderer','Painters');
print('-depsc','-tiff','-r300', '-painters', fullfile(pwd, 'example_gait_zscore.eps'));

