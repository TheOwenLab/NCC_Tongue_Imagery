% Code to perform analyse GLM stats
clear

% load GLM stats: 1st level analysis
%load GLM_processed_patient_MI_005_05_detrend_Aug27.mat
load GLM_processed_patient.mat

% Create Vector of Contrast:
% The vector of contrast is huge to mimic the case in which
% all available regressors were good.
% But only the first entries are important.

C = zeros(1,17);

C(1) = 1;


% List with all Short-channels
SSlist = [2 4 6 8 10 12 14 16];

% Get Stats from the First Level
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
    
    Act{cnt_sub} = find(pFDRHbO<0.05 & pFDRHbR<0.05 & ...
        T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);
    
    % Save activated channels following the original index
    Act_2{Nsub} = find(pFDRHbO<0.05 & pFDRHbR<0.05 & ...
        T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);

    % isolate deactivators
    Act_3{Nsub} = find((pFDRHbO<0.05 & pFDRHbR<0.05 & ...
        T(cnt_sub,:,1)<0 & T(cnt_sub,:,2)>0));
    
    % either activator or deactivator
    Act_4{Nsub} = sort([Act_2{Nsub},Act_3{Nsub}]);
    
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

%% Writing the T values for each participant to an excel sheet in matlab 

columnNamePHbO = []; columnNamePHbR = [];
for roi = 1:nROI
    columnNamePHbO = [columnNamePHbO,'Channel_' + string(roi) + '_' + 'HbOp'];
    columnNamePHbR = [columnNamePHbR,'Channel_' + string(roi) + '_' + 'HbRp'];
end

columnNameTHbO = []; columnNameTHbR = [];
for roi = 1:nROI
    columnNameTHbO = [columnNameTHbO,'Channel_' + string(roi) + '_' + 'HbOt'];
    columnNameTHbR = [columnNameTHbR,'Channel_' + string(roi) + '_' + 'HbRt'];
end

% also store patient names
load('patient_index.mat');
clean_subs = subnames(AvailableParticipants);

results = array2table([AvailableParticipants',P(:,:,1),P(:,:,2),T(:,ROI_group,1),T(:,ROI_group,2)],...
    'VariableNames',['Participants', columnNamePHbO,columnNamePHbR,columnNameTHbO,columnNameTHbR]);

results{:,'ID'} = clean_subs';
%% check the proportion of patients with activated channels
% any channel
Act_clean = Act_2(AvailableParticipants);
fprintf('Proportion of patients with any channels with tongue imagery %.4f \n', sum(~cellfun(@isempty,Act_clean))/length(Act_clean))
% which participants
disp('Which participants') 
SigParts = find(~cellfun(@isempty,Act_clean))

FrontalROI = [5,7]; MotorROI = [1,3,9,11]; ParietalROI = [13,15];
ClassROI = sort([MotorROI,FrontalROI]);

% not just any channel  
Act_Act = Act_clean(find(~cellfun(@isempty,Act_clean)));

MotorMask = zeros(length(Act_Act),1);
% which participants show non-parietal activation
for sp = 1:length(MotorMask)
    inter = intersect(Act_Act{sp},ClassROI);
    MotorMask(sp) = ~isempty(inter);
end

fprintf('Proportion of patients with any channels with tongue imagery in frontal or motor roi %.4f \n', sum(MotorMask)/length(Act_clean))

disp('Which Participants')
SigParts(MotorMask==1)


% and what about the deactivators?
Act_clean_deact = Act_3(AvailableParticipants);
fprintf('Proportion of patients with any channels with tongue imagery deactivation %.4f \n', ....
    sum(~cellfun(@isempty,Act_clean_deact))/length(Act_clean_deact))
% which participants
disp('Which participants') 
SigPartsDeact = find(~cellfun(@isempty,Act_clean_deact));

% not just any channel  
Act_Deact = Act_clean_deact(find(~cellfun(@isempty,Act_clean_deact)));

MotorMaskDeact = zeros(length(Act_Deact),1);
% which participants show non-parietal activation
for sp = 1:length(MotorMaskDeact)
    inter = intersect(Act_Deact{sp},ClassROI);
    MotorMaskDeact(sp) = ~isempty(inter);
end

fprintf('Proportion of patients with any channels with deactivation during tongue imagery in frontal or motor roi %.4f \n', sum(MotorMaskDeact)/length(Act_clean_deact))

disp('Which Participants')
SigPartsDeact(MotorMaskDeact==1)

% total proportion
disp('Proportion of patients who show either activation or deactivation')
length(unique(sort([SigParts(MotorMask==1);SigPartsDeact(MotorMaskDeact==1)])))/length(AvailableParticipants)



% store these results in the results excel
AnySig = zeros(length(AvailableParticipants),1); AnySig(SigParts) = 1;
NotParSig = zeros(length(AvailableParticipants),1); NotParSig(SigParts(MotorMask==1)) = 1;
AnySigDeact = zeros(length(AvailableParticipants),1); AnySigDeact(SigPartsDeact) = 1;
NotParSigDeact = zeros(length(AvailableParticipants),1); NotParSigDeact(SigPartsDeact(MotorMaskDeact==1)) = 1;

% add to results table
results{:,'AnySig'} = AnySig;
results{:,'NotParSig'} = NotParSig;
results{:,'AnySigDeact'} = AnySigDeact;
results{:,'NotParSigDeact'} = NotParSigDeact;

writetable(results,'resultsfdr.xlsx')
%% plot results

Act_within_Roi = cell(nPart,1); Act_within_Roi_HbO = cell(nPart,1);
Act_within_Roi_HbR = cell(nPart,1); 


for Nsub = AvailableParticipants

    % Activated channels within the ROI
    Act_within_Roi{Nsub} = ...
        intersect(Act_2{Nsub},ROI_group);


    Act_within_Roi_HbO{Nsub} = ...
        intersect(Act_HbO{Nsub},ROI_group);

    Act_within_Roi_HbR{Nsub} = ...
        intersect(Act_HbR{Nsub},ROI_group);

    % plot Activated Channels
    Imap = zeros(16,1);
    Imap(Act_within_Roi{Nsub}) = 1;

    lambda = 0.000005;
     Pardis_Create_3D_Plot_projection_MC_Androu_2...
         (Imap,lambda,[-1 1],[],jet(256));

    print(['Patient_' num2str(Nsub) '_Top_Contrast_'],'-dpng','-r300');
    close
        set(gcf,'color','w');
        set(gcf, 'InvertHardcopy', 'off');

    print(['Patient_' num2str(Nsub) '_Left_Contrast_'],'-dpng','-r300');
    close 

    %print(['Patient_' num2str(Nsub) '_Right_Contrast_'],'-dpng','-r300');
    close all;

end


