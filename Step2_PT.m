% Code to perform first level GLM analysis for patient data
clear

% load processed data
% Set local folder path
%load('processed_patient_data_individual_dpf_MI_Aug27');
%load('processed_patient_data_individual_dpf_MI_July18')
load('processed_patient_data');


clear gdata;

% Short Channels List
SSlist = [2 4 6 8 10 12 14 16];
GoodSC_Patients_MI;

% GLM option and Stim duration
opt_GLM = 1;

duration{1} = 15;

nPart = length(data);
start = 3; 

% exclude short channels from bad channels list as the definition of
% whether the channel is bad comes from whether a visible heart rate is
% present.
for part = 1:nPart
    tempchan =  BadChan{part};
    [r,c] = find(tempchan == SSlist);
    BadChan{part}(r) = [];

end

AvailableParticipants = start:nPart;
% drop patients without a complete dataset
AvailableParticipants(AvailableParticipants == 9) = [];
AvailableParticipants(AvailableParticipants == 20) = [];
AvailableParticipants(AvailableParticipants == 33) = [];
AvailableParticipants(AvailableParticipants == 52) = [];

AvailableParticipants = AvailableParticipants(1:end);

% Run analysis for all participants
% just used for skipping patient one for now since they have no good short channels and this may be causing an error

% preallocate variables since this code would not work if last participant
% is cut from the data

data_filtered_sc = cell(nPart,1); beta_R = cell(nPart,1);
covb_R = cell(nPart,1); p_R = cell(nPart,1);

for Nsub = AvailableParticipants
    
    r =  data{Nsub}{1};
    
    % Low-pass filter the data
    r.dc = r.BPFilter(r.dc, [0 0.2]);
    
    % detrend the data
    %r.dc(:,:,1) = detrend(r.dc(:,:,1));
    %r.dc(:,:,2) = detrend(r.dc(:,:,2));
    %r.dc(:,:,3) = detrend(r.dc(:,:,3));
    
    % Remove data before the first "Rest" trigger (Not done for patients
    % since triggers start at t=0)
    
    % remove task between the blocks
    [stimrow, stimcol] = find(r.s);

    restrow = stimrow(stimcol==1); restcol = stimcol(stimcol==1);
    trow = stimrow(stimcol==2); tcol = stimcol(stimcol==2);

    % first drop time before task starts (instruction period)
    start_drop = 1:restrow(1)-1; % perhaps cut first rest period since it isn't modelled by hrf?

    % then drop first break
    b1_start = ceil(trow(8) + duration{1} * r.SD.f);
    b1_drop = b1_start:restrow(9)-1;

    % second break
    b2_start = ceil(trow(16) + duration{1} * r.SD.f);
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
    
    % Perform SC regression like RS data
    r_filtered  = r;
    
    [r_filtered.dc,Stats] = ...
        r_filtered.PerformPhysiologyRegression...
        (GoodSC{Nsub},[],100,opt_GLM,...
        0,0);
    
    data_filtered_sc{Nsub}{1} = r_filtered;
    
    clear r_filtered;
    
    % Remove bad channels from hemoglobin time series
    % before perfoming the GLM analysis to infer activated 
    %channels
    r.dc(:,BadChan{Nsub},:) = nan;
        
    
    [beta_R{Nsub},p_R{Nsub},covb_R{Nsub}] = r.GLM_Prewhitening...
        (duration,GoodSC{Nsub},[],[],[],[],opt_GLM);
    %GLM_Prewhitening
    % Update data
    data{Nsub}{1} = r;
    
    
    %figure
    %plot(stim)
    %hold on 
    %plot(squeeze(r.dc(:,1,1)))
    %hold off
    
end

save('GLM_processed_patient','beta_R','p_R',...
    'covb_R','BadChan','data',...
    'data_filtered_sc',...
    'AvailableParticipants','SSlist');







