%%% Step 1: Analyze MI Data Patients
%%%
%%% In this first step, only the Hemoglobin time-series
%%% are computed
clear

% Available Participants:
% We are just including participants that have at least 
% one good short channel
AvailableParticipants = 1:52; 

% load data (Controls) -
% Set local folder path
load('PatientData.mat');

% Short Channels List
SSlist = [2 4 6 8 10 12 14 16]; % Updated to fit Pardis's short channels 

% Spline values
IndividualSplineValuesPatients_MI;

% Find Channels that do not have acceptable SNR
preproc.SNR_threshold = 8;

% Run for Quality Check for all Subjects
for Nsub = AvailableParticipants
    %%% ******* Step 1: NIRS DATA *********
    
    % Take NIRS data from specific Subject and Run
    obj = gdata{1,Nsub};
    
    % find bad channels
    Baseline = mean(obj.d);
    obj.d = detrend(obj.d)+Baseline;
    obj = obj.MarkBadChannels(preproc);
    BadChan{Nsub} = obj.SD.BadChannels;
    
    clear obj Baseline;
end

% Calculate proper dpf for each participant by taking into
% account their ages
Age_Patients_MI;

% extratc Lambda for each participant
for Nsub = AvailableParticipants    
    
   Wavelegnt_Sub{Nsub} =  gdata{1,Nsub}.SD.Lambda; 
   
end

alpha=223.3; 
beta=0.05624;
gamma=0.8493; 
delta=-5.723*10^-7;
epsilon=0.001245; 
zeta=-0.9025; 

for Nsub = AvailableParticipants    
    
    cnt=0;
    for wave_number= Wavelegnt_Sub{Nsub}
        cnt = cnt+1;
        dpf{Nsub}(cnt)=...
            alpha+beta*(Age(Nsub)^gamma)+delta*(wave_number^3)+...
         epsilon*(wave_number^2)+zeta*wave_number; 
     
    end
    
end


%Run preprocessing for all participants
for Nsub = AvailableParticipants
    
    % Get cw_nirs object
    r =  gdata{1,Nsub};
    
    % Compute Optical Density
    dOD = r.Convert2OD();
    
    r.SD.MeasListAct = ones(size(dOD,2),1);

    Opt_TDDR = 1;
    dodTDDR = hmrMotionCorrectTDDR_adapted(dOD,r.SD,r.SD.f,Opt_TDDR);
    
    % Estimate concentration
    r.dc = r.Convert2Conc(dodTDDR,dpf{Nsub});
    
    % Save data for further analysis
    data{Nsub}{1} = r;
    
    clear r;
    
end

% Save Hemoglobin time series for further analysis
save('processed_patient_data',...
    'data','BadChan','AvailableParticipants','gdata');

