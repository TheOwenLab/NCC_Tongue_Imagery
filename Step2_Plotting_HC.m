


% Code to perform first level GLM analysis for patient data
clear


load('preprocessed_data');


clear gdata;

% Short Channels List
SSlist = [2 4 6 8 10 12 14 16];
GoodSC_Controls_MI;

% GLM option and Stim duration
opt_GLM = 1;

restDur = 14; taskDur = 15;

% determine which channels to plot
plotROI = 1:2:16; MotorROI = [1,3,9,11]; MotorROIPlot = [1,2,5,6];

duration{1} = taskDur;

nPart = length(data);
start = 1; 

% exclude short channels from bad channels list as the definition of
% whether the channel is bad comes from whether a visible heart rate is
% present.
for part = 1:nPart
    tempchan =  BadChan{part};
    [r,c] = find(tempchan == SSlist);
    BadChan{part}(r) = [];

end

AvailableParticipants = start:nPart;

% load in subject (folder) names
load control_index.mat
clean_subs = subnames;

% Run analysis for all participants

% preallocate variables since this code would not work if last participant
% is cut from the data

data_filtered_sc = cell(nPart,1); beta_R = cell(nPart,1);
covb_R = cell(nPart,1); p_R = cell(nPart,1);
stim_vector_group = cell(nPart,1);

GroupRest = cell(nPart,1); GroupMI = cell(nPart,1);

% adjust this value to set lim for individual and group plots
 lim = [-0.35,0.35]; %individual plots
 group_lim = [-0.2,0.2]; % group limit

for Nsub = AvailableParticipants
    r =  data{Nsub}{1};

    % Low-pass filter the data
    r.dc = r.BPFilter(r.dc, [0 0.2]);
    
    % detrend the data
    %r.dc(:,:,1) = detrend(r.dc(:,:,1));
    %r.dc(:,:,2) = detrend(r.dc(:,:,2));
    %r.dc(:,:,3) = detrend(r.dc(:,:,3));
    
    % Remove data before the first "Rest" trigger 
    
    % remove task between the blocks
    [stimrow, stimcol] = find(r.s);

    restrow = stimrow(stimcol==1); restcol = stimcol(stimcol==1);
    trow = stimrow(stimcol==2); tcol = stimcol(stimcol==2);

    % first drop time before task starts (instruction period)
    start_drop = 1:restrow(1)-1; 

    % then drop first break
    b1_start = round(trow(8) + duration{1} * r.SD.f);
    b1_drop = b1_start:restrow(9)-1;

    % second break
    b2_start = round(trow(16) + duration{1} * r.SD.f);
    b2_drop = b2_start:restrow(17)-1;

    % drop end
    end_drop = ceil(trow(end) + duration{1} * r.SD.f);
    end_drop = end_drop:length(r.s);

    % readjust data to remove datapoints outside of task
    drop_idx = [start_drop, b1_drop, b2_drop, end_drop];

    r.s(drop_idx,:) = [];
    r.t(drop_idx,:) = [];
    r.d(drop_idx,:) = [];
    r.dc(drop_idx,:,:) = [];
    
    % readjust timing variable
    r.t = 1/r.SD.f:1/r.SD.f:(1/r.SD.f)*length(r.t);
    r.t = r.t';

    % Remove the first type of trigger
    stim = r.s; % for plotting later
    r.s(:,1) = [];
    
    % adjust the stim vector to one that is better for plotting
    % find conditions
    mi_start = find(r.s);
    stimVector = zeros(length(r.s),1);
    
    % for each condition
    for trial = 1:length(mi_start)
        stimVector(mi_start(trial):mi_start(trial)+floor(r.SD.f * duration{1})) = 1;
    end
    
    stim_vector_group = stimVector;

       % Perform SC regression like RS data
    r_filtered  = r;
    
    [r_filtered.dc, Stats] = ...
        r_filtered.PerformPhysiologyRegression...
        (GoodSC{Nsub}, [], 100, 4, ...
        0, 0);
    
    data_filtered_sc{Nsub}{1} = r_filtered;

          % Initialize variables
    epoch_data = data_filtered_sc{Nsub}{1}.dc(:,plotROI,:);
    [total_samples, nROI, nHb] = size(epoch_data);
    rest_samples = round(restDur * r.SD.f); % Number of samples for the rest condition
    mi_samples = round(taskDur * r.SD.f); % Number of samples for the motor imagery condition
   
    % Find the indices where the conditions change
    condition_changes = find(diff(stimVector) ~= 0) + 1;
    condition_changes = [1; condition_changes]; % Add the first trial
    num_trials = length(condition_changes) / 2; % Each trial has both rest and mi
   
    % Preallocate matrices for rest and motor imagery conditions
    RestEpochs = zeros(num_trials, rest_samples, nROI, nHb);
    MiEpochs = zeros(num_trials, mi_samples, nROI, nHb);
   
    % Loop through each trial and extract the data based on stimVector
    for i = 1:num_trials
        rest_start = condition_changes((i-1)*2 + 1);
        mi_start = condition_changes((i-1)*2 + 2);
       
        RestEpochs(i, :, :, :) = epoch_data(rest_start:rest_start+rest_samples-1, :, :);
        MiEpochs(i, :, :, :) = epoch_data(mi_start:mi_start+mi_samples-1, :, :);
    end

    % Save Epochs
    GroupRest{Nsub} = RestEpochs;
    GroupMI{Nsub} = MiEpochs;

    % Get the sampling frequency
    fs = r.SD.f;
   
    % Convert sample space to time (seconds)
    time_rest = (0:rest_samples-1) / fs;
    time_mi = (0:mi_samples-1) / fs;
    
		% calculate limits across plots
	%max_lim = max(max(abs(mean(MiEpochs(:, :, :, :),1)),[],'all'),...
      %  max(abs(mean(RestEpochs(:, :, :, :),1)),[],'all'));

   % max_error = max(max(std(MiEpochs(:, :, :, 1),[],1),[],'all'), ...
     %    max(std(RestEpochs(:, :, :, 1),[],1),[],'all'));

    %lim = [-max_lim - max_error/sqrt(num_trials), max_lim + max_error/sqrt(num_trials)];
   

    % Create a new figure for each set of plots
    pdfFileName = ['Participant_' clean_subs{Nsub} '_Plots.pdf'];
    all_figures = []; % Collect figure handles
   
    % Plot each channel in the ROI for both HbO and HbR, overlaying rest and motor imagery
    for roi = 1:nROI
        fig = figure('Visible', 'off'); % Create a new figure and keep it hidden for each ROI
        sgtitle(['Participant ' clean_subs{Nsub} ' - Channel ' num2str(plotROI(roi))]);
       
        % first subplot MI
        subplot(1, 2, 1);
        hold on
        colour_map = ['r','b'];
		
		% calculate limits across conditions but the same for each ROI 
		%lim = [min(mean(MiEpochs(:, :, roi, :),'all'),mean(RestEpochs(:, :, roi, :),'all')),...
		%max(mean(MiEpochs(:, :, roi, :),'all'),mean(RestEpochs(:, :, roi, :),'all'))];
		
        for hb = 1:nHb
            % Calculate mean and SEM for the motor imagery condition
            avg_mi = squeeze(mean(MiEpochs(:, :, roi, hb), 1));
            sem_mi = squeeze(std(MiEpochs(:, :, roi, hb), 0, 1)) / sqrt(num_trials);
                     
            % Plot motor imagery condition with error bars
            shadedErrorBar(time_mi, avg_mi, sem_mi, 'lineProps', {colour_map(hb), 'LineWidth', 1.5});            
        end
        title(['Tongue Imagery']);
        xlabel('Time (s)');
        ylim(lim);
        ylabel('Concentration Change (\DeltaHb)');
        legend('HbO', 'HbR');
        hold off
         
        % second subplot REST
        subplot(1, 2, 2);
        hold on

        for hb = 1:nHb
            % Calculate mean and SEM for the rest condition
            avg_rest = squeeze(mean(RestEpochs(:, :, roi, hb), 1));
            sem_rest = squeeze(std(RestEpochs(:, :, roi, hb), 0, 1)) / sqrt(num_trials);
            % Plot rest condition with error bars
            shadedErrorBar(time_rest, avg_rest, sem_rest, 'lineProps', {colour_map(hb), 'LineWidth', 1.5});
        end
       
        % Set condition labels
        title(['Rest']);
        xlabel('Time (s)');
        ylabel('Concentration Change (\DeltaHb)');
        ylim(lim);
        legend('HbO', 'HbR');
        hold off

   
        all_figures = [all_figures fig];
    end

       
    % Plot the average across all channels in the ROI
    fig = figure('Visible', 'off');
    sgtitle(['Participant ' clean_subs{Nsub}  ' - Average Across Channels']);
   
    % first plot motor imagery
    subplot(1, 2, 1);
    hold on
    for hb = 1:nHb
        avg_mi = squeeze(mean(mean(MiEpochs(:, :, :, hb), 1), 3));
        sem_mi = squeeze(std(mean(MiEpochs(:, :, :, hb), 3), 0, 1)) / sqrt(num_trials);
        shadedErrorBar(time_mi, avg_mi, sem_mi, 'lineProps', {colour_map(hb), 'LineWidth', 1.5});
    end
    title(['Tongue Imagery']);
    xlabel('Time (s)');
    ylabel('Concentration Change (\DeltaHb)');
    ylim(lim);
    legend('HbO', 'HbR');
    % then rest
    subplot(1, 2, 2);
    hold on

    for hb = 1:nHb
        avg_rest = squeeze(mean(mean(RestEpochs(:, :, :, hb), 1), 3));
        sem_rest = squeeze(std(mean(RestEpochs(:, :, :, hb), 3), 0, 1)) / sqrt(num_trials);
        shadedErrorBar(time_rest, avg_rest, sem_rest, 'lineProps', {colour_map(hb), 'LineWidth', 1.5});
    end
   
    title(['Rest']);
    xlabel('Time (s)');
    ylabel('Concentration Change (\DeltaHb)');
    ylim(lim);
    legend('HbO', 'HbR');
   
    all_figures = [all_figures fig];
   
    % Plot the average across channels in the MotorROIPlot
    fig = figure('Visible', 'off');
    sgtitle(['Participant ' clean_subs{Nsub}  ' - Average Across Motor ROI Channels']);
   
    motor_roi_indices = MotorROIPlot;
   
     % first plot motor imagery
    subplot(1, 2, 1);
    hold on
    for hb = 1:nHb
        avg_mi = squeeze(mean(mean(MiEpochs(:, :, motor_roi_indices, hb), 1), 3));
        sem_mi = squeeze(std(mean(MiEpochs(:, :, motor_roi_indices, hb), 3), 0, 1)) / sqrt(num_trials);
        shadedErrorBar(time_mi, avg_mi, sem_mi, 'lineProps', {colour_map(hb), 'LineWidth', 1.5});
    end
    title(['Tongue Imagery']);
    xlabel('Time (s)');
    ylabel('Concentration Change (\DeltaHb)');
    ylim(lim);
    legend('HbO', 'HbR');
    
	% then rest
    subplot(1, 2, 2);
    hold on
    for hb = 1:nHb
        avg_rest = squeeze(mean(mean(RestEpochs(:, :, motor_roi_indices, hb), 1), 3));
        sem_rest = squeeze(std(mean(RestEpochs(:, :, motor_roi_indices, hb), 3), 0, 1)) / sqrt(num_trials);
        shadedErrorBar(time_rest, avg_rest, sem_rest, 'lineProps', {colour_map(hb), 'LineWidth', 1.5});
    end
   
    title(['Rest']);
    xlabel('Time (s)');
    ylabel('Concentration Change (\DeltaHb)');
    ylim(lim);
    legend('HbO', 'HbR');

    all_figures = [all_figures fig];
   
    % Plot each individual trial for the ROI
    for roi = 1:nROI
        fig = figure('Visible', 'off');
        sgtitle(['Participant ' clean_subs{Nsub}  ' - Individual Trials for ROI ' num2str(plotROI(roi))]);

        for hb = 1:nHb
            subplot(1, nHb, hb);
            % Plot each rest trial for the current ROI and hemoglobin type
            for trial = 1:num_trials
                plot(time_rest, squeeze(RestEpochs(trial, :, roi, hb)), 'Color', [0, 0, 1, 0.3]); % Blue with transparency for rest
                hold on;
            end
           
            % Plot each motor imagery trial for the current ROI and hemoglobin type
            for trial = 1:num_trials
                plot(time_mi, squeeze(MiEpochs(trial, :, roi, hb)), 'Color', [1, 0, 0, 0.3]); % Red with transparency for motor imagery
                hold on;
            end
           
            hb_type = {'HbO', 'HbR'};
            title([hb_type{hb} ' (ROI ' num2str(plotROI(roi)) ')']);
            xlabel('Time (s)');
            ylabel('Concentration Change (\DeltaHb)');
            legend({'Rest Trials', 'Motor Imagery Trials'}, 'Location', 'best');
        end
        all_figures = [all_figures fig];
    end

    % Save all figures to a single PDF
    for i = 1:length(all_figures)
        exportgraphics(all_figures(i), pdfFileName, 'Append', true);
    end
   
    % Close all figures after saving
    close(all_figures);
            
    
    [beta_R{Nsub},p_R{Nsub},covb_R{Nsub}] = r.GLM_Prewhitening...
        (duration,GoodSC{Nsub},[],[],[],[],opt_GLM);
    
	%GLM_Prewhitening
    % Update data
    data{Nsub}{1} = r;
    
end

%% plot group level results
% Get the number of participants, ROIs, and hemoglobin types
nPart = length(GroupRest);
[nTrials, rest_samples, nROI, nHb] = size(GroupRest{1});

% Get the sampling frequency (assuming consistent across participants)
fs = r.SD.f;

% Convert sample space to time (seconds)
time_rest = (0:rest_samples-1) / fs;
time_mi = (0:size(GroupMI{1}, 2)-1) / fs;  % Use GroupMI samples length for time

% Preallocate arrays for group-level data
group_avg_RestEpochs = cell(1, nROI);
group_avg_MiEpochs = cell(1, nROI);

for roi = 1:nROI
    % Preallocate arrays for group-level averaging
    group_RestData = zeros(nPart, rest_samples, nHb);
    group_MiData = zeros(nPart, length(time_mi), nHb);
    
    % Collect and average data across participants for each ROI
    for Nsub = 1:nPart
        group_RestData(Nsub, :, :) = squeeze(mean(GroupRest{Nsub}(:, :, roi, :), 1));
        group_MiData(Nsub, :, :) = squeeze(mean(GroupMI{Nsub}(:, :, roi, :), 1));
    end
    
    % Calculate the group mean across participants for each ROI
    group_avg_RestEpochs{roi} = squeeze(mean(group_RestData, 1));
    group_avg_MiEpochs{roi} = squeeze(mean(group_MiData, 1));
end

% Create PDF to save group plots
group_pdfFileName = 'Group_Level_Plots.pdf';
all_group_figures = []; % Collect figure handles
% Group-level average across all channels
fig = figure('Visible', 'off');
sgtitle('Group Level - Average Across All Channels');

% Calculate the group-level mean and SEM across all ROIs
avg_group_rest = mean(cat(3, group_avg_RestEpochs{:}), 3);
sem_group_rest = std(cat(3, group_avg_RestEpochs{:}), 0, 3) / sqrt(nROI);

avg_group_mi = mean(cat(3, group_avg_MiEpochs{:}), 3);
sem_group_mi = std(cat(3, group_avg_MiEpochs{:}), 0, 3) / sqrt(nROI);

% First panel: HbO and HbR activity for Motor Imagery condition
subplot(1, 2, 1);
hold on;
shadedErrorBar(time_mi, avg_group_mi(:, 1), sem_group_mi(:, 1), 'lineProps', {'r', 'LineWidth', 1.5}); % HbO
shadedErrorBar(time_mi, avg_group_mi(:, 2), sem_group_mi(:, 2), 'lineProps', {'b', 'LineWidth', 1.5}); % HbR
title('Motor Imagery (All Channels)');
xlabel('Time (s)');
ylim(group_lim)
ylabel('Concentration Change (\DeltaHb)');
legend('HbO', 'HbR');
hold off;

% Second panel: HbO and HbR activity for Rest condition
subplot(1, 2, 2);
hold on;
shadedErrorBar(time_rest, avg_group_rest(:, 1), sem_group_rest(:, 1), 'lineProps', {'r', 'LineWidth', 1.5}); % HbO
shadedErrorBar(time_rest, avg_group_rest(:, 2), sem_group_rest(:, 2), 'lineProps', {'b', 'LineWidth', 1.5}); % HbR
title('Rest (All Channels)');
xlabel('Time (s)');
ylim(group_lim)
ylabel('Concentration Change (\DeltaHb)');
legend('HbO', 'HbR');
hold off;

all_group_figures = [all_group_figures fig];

% Group-level average across channels within Motor ROI
fig = figure('Visible', 'off');
sgtitle('Group Level - Average Across Motor ROI Channels');

% Indices of Motor ROI
motor_roi_indices = MotorROIPlot;

% Calculate group-level mean and SEM for Motor ROI
group_MotorRest = zeros(nPart, rest_samples, nHb);
group_MotorMi = zeros(nPart, length(time_mi), nHb);

for Nsub = 1:nPart
    group_MotorRest(Nsub, :, :) = squeeze(mean(mean(GroupRest{Nsub}(:, :, motor_roi_indices, :), 1), 3));
    group_MotorMi(Nsub, :, :) = squeeze(mean(mean(GroupMI{Nsub}(:, :, motor_roi_indices, :), 1), 3));
end

avg_group_MotorRest = squeeze(mean(group_MotorRest, 1));
sem_group_MotorRest = squeeze(std(group_MotorRest, 0, 1)) / sqrt(nPart);

avg_group_MotorMi = squeeze(mean(group_MotorMi, 1));
sem_group_MotorMi = squeeze(std(group_MotorMi, 0, 1)) / sqrt(nPart);

% First panel: HbO and HbR activity for Motor Imagery condition
subplot(1, 2, 1);
hold on;
shadedErrorBar(time_mi, avg_group_MotorMi(:, 1), sem_group_MotorMi(:, 1), 'lineProps', {'r', 'LineWidth', 1.5}); % HbO
shadedErrorBar(time_mi, avg_group_MotorMi(:, 2), sem_group_MotorMi(:, 2), 'lineProps', {'b', 'LineWidth', 1.5}); % HbR
title('Motor Imagery (Motor ROI)');
xlabel('Time (s)');
ylim(group_lim)
ylabel('Concentration Change (\DeltaHb)');
legend('HbO', 'HbR');
hold off;

% Second panel: HbO and HbR activity for Rest condition
subplot(1, 2, 2);
hold on;
shadedErrorBar(time_rest, avg_group_MotorRest(:, 1), sem_group_MotorRest(:, 1), 'lineProps', {'r', 'LineWidth', 1.5}); % HbO
shadedErrorBar(time_rest, avg_group_MotorRest(:, 2), sem_group_MotorRest(:, 2), 'lineProps', {'b', 'LineWidth', 1.5}); % HbR
title('Rest (Motor ROI)');
xlabel('Time (s)');
ylim(group_lim)
ylabel('Concentration Change (\DeltaHb)');
legend('HbO', 'HbR');
hold off;

all_group_figures = [all_group_figures fig];

% Save all group-level figures to a single PDF
for i = 1:length(all_group_figures)
    exportgraphics(all_group_figures(i), group_pdfFileName, 'Append', true);
end

% Close all figures after saving
close(all_group_figures);