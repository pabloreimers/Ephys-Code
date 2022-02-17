
ephysSettings

basePath    = settings.mainDataDir;
date        = input('Date: ','s');
genotype    = input('Genotype: ','s');
nfly        = input('Fly ID: ','s');
cell_num    = input('Cell ID: ','s');
exp_name    = input('Exp Name: ','s');

trial_repeats   = input('trial repeats:  ');
start_length    = input('start length (s):   ');
stim_length     = input('stim length (s):    ');
end_length      = input('end length (s):   ');
inR_log         = input('calculate input resistance? (logical):  ');

[rawData, trialData, trialMeta] = acquireSimpleTrial(trial_repeats, start_length, stim_length, end_length, inR_log);

if input('save params? (log): ')
    allData.rawData   = rawData;
    allData.trialData = trialData;
    allData.trialMeta = trialMeta;
    foldername  = sprintf('%s\\%s\\%s\\Fly %s\\Cell %s\\%s\\',settings.mainDataDir,date,genotype,nfly,cell_num,exp_name);
    num_file   = length(dir(foldername)) - 1;
    save_dir    = [foldername,num2str(num_file),'\'];
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    save([save_dir,'\allData.mat'],'allData');
    fprintf('\n********** Done **********\n')
end