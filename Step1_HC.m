%%% Step 1: Convert and Preprocess Data
%%%

clear

% load data (Controls) -
% Set local folder path
load('DataRaw_HC.mat');

% adjust the final number to be the total number of participants (with completed datasets)
% assumed right now to be the total number of participants in DataRaw_HC
AvailableParticipants = 1:length(gdata);  

% Short Channels List
SSlist = [2 4 6 8 10 12 14 16]; % Updated to fit Pardis's short channels 

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
Age_Controls_MI;

% extract Lambda for each participant
for Nsub = AvailableParticipants    
    
   Wavelength_Sub{Nsub} =  gdata{1,Nsub}.SD.Lambda; 
   
end

% parameters for computing the DPF
alpha=223.3; 
beta=0.05624;
gamma=0.8493; 
delta=-5.723*10^-7;
epsilon=0.001245; 
zeta=-0.9025; 

for Nsub = AvailableParticipants    
    
    cnt=0;
    for wave_number= Wavelength_Sub{Nsub}
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
    
    % apply TDDR motion correction
    Opt_TDDR = 1;
    dodTDDR = hmrMotionCorrectTDDR_adapted(dOD,r.SD,r.SD.f,Opt_TDDR);
    
    % Estimate concentration
    r.dc = r.Convert2Conc(dodTDDR,dpf{Nsub});
    
    % Save data for further analysis
    data{Nsub}{1} = r;
    
    clear r;
    
end

% Save Hemoglobin time series for further analysis
save('preprocessed_data',...
    'data','BadChan','AvailableParticipants','gdata');

