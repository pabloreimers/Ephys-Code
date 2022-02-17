
ephysSettings

basePath    = settings.mainDataDir;
date        = input('Date: ','s');
nfly        = input('Fly ID: ','s');
cell        = input('Cell ID: ','s');

trial_repeats   = input('trial repeats:  ');
trial_length    = input('trial length (s):   ');
inR_log         = input('calculate input resistance? (logical):  ');

[rawData, trialData, trialMeta] = acquireSimpleTrial(trial_repeats, trial_length, inR_log);

if input('save params? (log): ')
    allData.rawData   = rawData;
    allData.trialData = trialData;
    allData.trialMeta = trialMeta;
    foldername  = sprintf('%s\\%s\\Fly %s\\Cell %s\\',settings.mainDataDir,date,nfly,cell);
    num_file   = length(dir(foldername)) + 1;
    save_dir    = [foldername,'\',num2str(num_file),'\'];
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    save([save_dir,'\allData.mat'],'allData');
    fprintf('\n********** Done **********\n')
end