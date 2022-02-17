ephysSettings

basePath    = settings.mainDataDir;
date        = input('Date: ','s');
genotype    = input('Genotype: ','s');
nfly        = input('Fly ID: ','s');
cell_num    = input('Cell ID: ','s');
exp_name    = input('Exp Name: ','s');

%acquire meta data recordings. display where pipette should be, and key
%press to record the trial. if recording is good, ask user and proceed.
%otherwise, just redo.
% meta_flag = false;
% tmp       = false;
% while ~meta_flag
%     while ~tmp ; fprintf('bath?') ;      pause ; [~, meta.bath, ~]   = acquireSimpleTrial(1,.5,0,0,0); tmp = input('good?'); end; tmp = false;
%     while ~tmp ; fprintf('seal?') ;      pause ; [~, meta.seal, ~]   = acquireSimpleTrial(1,.5,0,0,0); tmp = input('good?'); end; tmp = false;
%     while ~tmp ; fprintf('patch?');      pause ; [~, meta.patch, ~]  = acquireSimpleTrial(1,.5,0,0,0); tmp = input('good?'); end; tmp = false;
%     while ~tmp ; fprintf('baseline?');   pause ; [~, meta.baseline, ~] = acquireSimpleTrial(1,10,0,0,0); tmp = input('good?'); end; tmp = false;
%     meta_flag = input('all good?  ');
% end
% 
% allData.meta    = meta;
trial_repeats   = input('trial repeats:  ');
start_length    = input('start length (s):   ');
stim_length     = input('stim length (s):    ');
end_length      = input('end length (s):   ');

disp('I-clamp')
pause
[~, baseline, ~]                        = acquireSimpleTrial(1,10,0,0,0);
pause
[rawData, trialData, trialMeta]         = acquireRunningTrial(trial_repeats, start_length, stim_length, end_length, true);
trialMeta.pipetteR                      = input('pipette resistance (MOhms): ');
trialMeta.sealR                         = input('seal resistance (GOhms): ');
trialMeta.access                        = input('access resistance (MOhms): ');
trialMeta.baseline                      = baseline;

if input('save params? (log): ')
    clear allData
    allData.rawData   = rawData;
    allData.trialData = trialData;
    allData.trialMeta = trialMeta;
    foldername  = sprintf('%s\\%s\\%s\\Fly %s\\Cell %s\\%s\\',settings.mainDataDir,date,genotype,nfly,cell_num,exp_name);
    num_file   = length(dir(foldername)) - 1;
    save_dir    = [foldername,num2str(num_file),'\'];
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    save([save_dir,'\allData.mat'],'allData','-v7.3');
    fprintf('\n********** Done **********\n')
end