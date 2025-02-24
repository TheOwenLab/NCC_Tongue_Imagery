% Code to perform first level GLM analysis for patient data
clear

% load processed data
load('preprocessed_data');


clear gdata;

% Short Channels List
SSlist = [2 4 6 8 10 12 14 16];
Chans = 1:16;
GoodSC_Controls_MI;

% GLM option and Stim duration
restDur = 14; taskDur = 15;
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

% Run analysis for all participants

% preallocate variables since this code would not work if last participant
% is cut from the data

% Define grid of parameter values
lowpass_values = [0, 0.005, 0.01];
highpass_values = [0.1, 0.2, 0.3, 0.4, 0.5, 1];
detrend_values = [false, true];
optGLM_values = [1, 2, 4];

% Generate all combinations of parameters
[lowpass_grid, highpass_grid, detrend_grid, optGLM_grid] = ndgrid(lowpass_values, highpass_values, detrend_values, optGLM_values);

% Preallocate result variables
num_combinations = numel(lowpass_grid);

Results.Act = cell(num_combinations,1); Results.Act_Prop = cell(num_combinations,1); 
Results.SigParts = cell(num_combinations,1); Results.Act_Prop_NonParietal = cell(num_combinations,1); 
Results.SigPartsNonParietal = cell(num_combinations,1); 
Results.Act_Prop_Deact = cell(num_combinations,1); 
Results.SigPartsDeact = cell(num_combinations,1); 

Results.Params = cell(num_combinations,1); 


data_filtered_sc = cell(nPart,1); beta_R = cell(nPart,1);
covb_R = cell(nPart,1); p_R = cell(nPart,1);
stim_vector_group = cell(nPart,1);

for hyper_comb = 1:num_combinations
    % Get parameter values for current combination
    lowpass = lowpass_grid(hyper_comb);
    highpass = highpass_grid(hyper_comb);
    detrend_param = detrend_grid(hyper_comb);
    optGLM = optGLM_grid(hyper_comb);

    for Nsub = AvailableParticipants
        r =  data{Nsub}{1};
    
        % Low-pass filter the data
        r.dc = r.BPFilter(r.dc, [lowpass highpass]);
        
        if detrend_param
            % detrend the data
            r.dc(:,:,1) = detrend(r.dc(:,:,1));
            r.dc(:,:,2) = detrend(r.dc(:,:,2));
            r.dc(:,:,3) = detrend(r.dc(:,:,3));
        end
        
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
    
        r.s(:,1) = [];
    
    
           % Perform SC regression like RS data
        r_filtered  = r;
        
        [r_filtered.dc, Stats] = ...
            r_filtered.PerformPhysiologyRegression...
            (GoodSC{Nsub}, [], 100, 4, ...
            0, 0);
        
        data_filtered_sc{Nsub}{1} = r_filtered;
        
        % Remove bad channels from hemoglobin time series
        % before perfoming the GLM analysis to infer activated 
        %channels
        r.dc(:,BadChan{Nsub},:) = nan;
            
        
        [beta_R{Nsub},p_R{Nsub},covb_R{Nsub}] = r.GLM_Prewhitening...
            (duration,GoodSC{Nsub},[],[],[],[],optGLM);
        
        
    end
    
    %% run second level with different preprocessing 
    
    % Create Vector of Contrast:
    C = zeros(1,17);
    C(1) = 1;
    
    % determine where the significance is coming from 
    FrontalROI = [5,7]; MotorROI = [1,3,9,11]; ParietalROI = [13,15];
    ClassROI = sort([MotorROI,FrontalROI]);
    
    % Get Stats from the First Level Analysis
    [B_k,Cov_k] = ...
        ExtractDataFromFirstLeveL...
        (beta_R,covb_R,C,SSlist,AvailableParticipants,...
        BadChan);
    
    % Compute T-Values for All Subjects
    T = B_k./sqrt(Cov_k);
    
    % Perform intra-subject analysis for the defined contrast.
    % This steps will evaluate which channels were activated for
    % each participants. This result will used to compute the
    % sensitivity results
    
    cnt_sub = 0;
    at_least_one = [];
    at_least_one_sub = [];
    
    % Load ROI group Results
    %name_Roi_to_load = ['ROI_group_health_MI'];
    
    %ROI_group = load(name_Roi_to_load);
    ROI_group = 1:2:16;
    nROI = length(ROI_group);
    P = nan(length(AvailableParticipants),length(ROI_group),2);
    
    % preallocate again due to bug if last participant is cut 
    nPart = length(beta_R);
    Act = cell(nPart,1); Act_2 = cell(nPart,1);
    Act_3 = cell(nPart,1); Act_4 = cell(nPart,1);
    Act_HbO = cell(nPart,1); Act_HbR = cell(nPart,1);
    
    
    for Nsub = AvailableParticipants
        
        % Counter
        cnt_sub = cnt_sub +1;
        
        % Get cw-nirs object to infer the degress of freedom
        r = data{Nsub}{1};
        
        % Compute degrees of freedom as the length of the time-series
        DegreeOfFreedom = size(r.dc,1) - 50;
        
        % Convert T to P
        % HbO
        p_valueHbO = 1-tcdf(abs(T(cnt_sub,:,1)),DegreeOfFreedom);
        
        % HbR
        p_valueHbR = 1-tcdf(abs(T(cnt_sub,:,2)),DegreeOfFreedom);
    
        % fdr correct
        [~,~,~,pFDRHbO] = fdr_bh(p_valueHbO);
        [~,~,~,pFDRHbR] = fdr_bh(p_valueHbR);
             
        % Save activated channels following the original index
        Act_2{Nsub} = find(pFDRHbO<0.05 & pFDRHbR<0.05 & ...
            T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);
    
        % isolate deactivators
        Act_3{Nsub} = find((pFDRHbO<0.05 & pFDRHbR<0.05 & ...
            T(cnt_sub,:,1)<0 & T(cnt_sub,:,2)>0));
        
        Act_HbO{Nsub} = find(pFDRHbO<0.05 & ...
            T(cnt_sub,:,1)>0);
        
        Act_HbR{Nsub} = find(pFDRHbR<0.05 & ...
            T(cnt_sub,:,2)<0);
         
        % store p-values
        P(cnt_sub,:,1) = pFDRHbO(ROI_group);
        P(cnt_sub,:,2) = pFDRHbR(ROI_group);
        
        
        %end of section
    
    
        % Check which participants had at least one activated channel
        if ~isempty(Act{cnt_sub})
            at_least_one = [at_least_one,cnt_sub];
            at_least_one_sub = [at_least_one_sub,Nsub];
        end
        
        clear p_valueHbO p_valueHbR R DegreeOfFreedom
    
    end

    % start to store results
    Results.Act{hyper_comb} = Act_2; 

    %% check the proportion of patients with activated channels
    % any channel
    Act_clean = Act_2(AvailableParticipants);
    Act_Prop = sum(~cellfun(@isempty,Act_clean))/length(Act_clean);
    Results.Act_Prop{hyper_comb} = Act_Prop;
   
    SigParts = find(~cellfun(@isempty,Act_clean));
    Results.SigParts{hyper_comb} = SigParts;
    
    % not just any channel  
    Act_Act = Act_clean(find(~cellfun(@isempty,Act_clean)));
    
    MotorMask = zeros(length(Act_Act),1);
    % which participants show non-parietal activation
    for sp = 1:length(MotorMask)
        inter = intersect(Act_Act{sp},ClassROI);
        MotorMask(sp) = ~isempty(inter);
    end
    
    Act_Prop_NonParietal = sum(MotorMask)/length(Act_clean);
    Results.Act_Prop_NonParietal{hyper_comb} = Act_Prop_NonParietal;
    Results.SigPartsNonParietal{hyper_comb} = SigParts(MotorMask==1);

    % and what about the deactivators?
    Act_clean_deact = Act_3(AvailableParticipants);
    Prop_Act_deact = sum(~cellfun(@isempty,Act_clean_deact))/length(Act_clean_deact);
    Results.Act_Prop_Deact{hyper_comb} = Prop_Act_deact;
    SigPartsDeact = find(~cellfun(@isempty,Act_clean_deact));
    Results.SigPartsDeact{hyper_comb} = SigPartsDeact;
    
    Results.Params{hyper_comb} = [lowpass,highpass,detrend_param,optGLM];

end

% save the results variable
figure
plot([Results.Act_Prop_NonParietal{:}])

% find the maximum parameters
prop_act = [Results.Act_Prop_NonParietal{:}];
max_prop = max(prop_act);
best_indices = find(prop_act == max_prop);

fprintf("Percentage of participants with significant activation in non parietal regions %7.3f\n",max_prop)

disp([Results.Params{best_indices}])

