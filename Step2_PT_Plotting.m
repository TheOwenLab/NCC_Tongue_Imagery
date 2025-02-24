


% Code to perform first level GLM analysis for patient data
clear


load('processed_patient_data');


clear gdata;

% Short Channels List
SSlist = [2 4 6 8 10 12 14 16];
GoodSC_Patients_MI;

% GLM option and Stim duration
opt_GLM = 1;

restDur = 14; taskDur = 15;
duration{1} = taskDur;


nPart = length(data);
start = 3; 


% adjust this value to set lim for individual and group plots
 lim = [-0.25,0.25]; %individual plots
 group_lim = [-0.2,0.2]; % group limit

% exclude short channels from bad channels list as the definition of
% whether the channel is bad comes from whether a visible heart rate is
% present.
for part = 1:nPart
    tempchan =  BadChan{part};
    [r,c] = find(tempchan == SSlist);
    BadChan{part}(r) = [];

end

% pateient 
AvailableParticipants = start:nPart;
AvailableParticipants(AvailableParticipants == 9) = [];
AvailableParticipants(AvailableParticipants == 20) = [];
AvailableParticipants(AvailableParticipants == 33) = [];
AvailableParticipants(AvailableParticipants == 52) = [];

restDur = 14; taskDur = 15;
plotROI = 1:2:16; MotorROI = [1,3,9,11]; MotorROIPlot = [1,2,5,6];

% load in participant ids
load('patient_index.mat'); % change for patient in patient case
%clean_subs = subnames(AvailableParticipants);
clean_subs = subnames;


% Run analysis for all participants

% preallocate variables since this code would not work if last participant
% is cut from the data

data_filtered_sc = cell(nPart,1); beta_R = cell(nPart,1);
covb_R = cell(nPart,1); p_R = cell(nPart,1);
stim_vector_group = cell(nPart,1);

GroupRest = cell(nPart,1); GroupMI = cell(nPart,1);

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







