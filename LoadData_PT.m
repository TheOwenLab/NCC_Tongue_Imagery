%% generate .lob file to use Mesquita pipeline

% generate participant filenames
filepath = 'Patient_Data'; % if the files are in a different folder

instrument = 'NIRScout (NIRx)';
restDur = 14;

% load in the folder names to keep track of the folder to ID correpondence
subfolders = dir(filepath);
subs = subfolders([subfolders.isdir]); %making sure only directories are pulled
subs = subs(~ismember({subs.name},{'.','..'})); % removing "." and ".." from subs
subnames = {subs.name};
nsubs = length(subnames);

% change directory to where the files are stored

cd(filepath); 
for ss=1:nsubs
    % enter ss's folder
    cd([subnames{ss}])
    disp("Working on Participant " + subnames{ss});
    
    % grab probe and hdrfile
    probefile = dir('*Info.mat');
    hdrfile = dir('*.hdr');
    
    % grab trigger file
    tgrfile = dir('*.tri');

    % load data (note function will only work for this montage)
    data = ReadNIRSAurora(instrument,hdrfile.name,probefile.name); 

    % loading function truncates sampling frequency
    data.SD.f = 10.172526041666666;

    %load triggers
    T = load(tgrfile.name,"-ASCII");

    % store group data in cell
    s = zeros(length(data.t),2);

    % only triggers for task on
    s(T(:,2),2) = 1;

    % infer triggers for rest as restDur seconds before task
    s(round(find(s(:,2)) - (restDur*data.SD.f)),1) = 1;
    data.s = s;
    
    % store data
    gdata{ss} = data;
    % move back to the folder containing the group subject data
    cd ..
end

cd ..

% if .lob or gdata is not appearing you need to add the original data files
% to path 
save('PatientData','gdata')
save("patient_index","subnames")