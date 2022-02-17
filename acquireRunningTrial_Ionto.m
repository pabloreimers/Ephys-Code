function [rawData, trialData, trialMeta] = acquireRunningTrial_Ionto(trial_repeats, start_length, stim_length, end_length, inR_log)
% Acquisition code for baseline recording, no stimulus provided

% INPUT
    % trialRepeats - number of trial repeats
    % trialLength - duration of each trial
    % inR - [0/1] calculate input resistance?
% OUTPUT 
    % trialData
    % trialMeta
    % rawData
    
%% intialize settings

daqreset;
ephysSettings;
DAQsettings;

%% Add output channels
% create vector of analog outputs, the output will be written at each sample time for stim or no stim(however many seconds per however many samples per second)
tot_length = start_length + stim_length + end_length;
output      = [zeros(ceil(tot_length * niIO.Rate),1)]; 

A0 = addoutput(niIO,settings.devID,'ao0','Voltage'); %add an output channel to our NI breakout board
A0.Name = 'Shutter';

A1 = addoutput(niIO,settings.devID,'ao1','Voltage');
A1.Name = 'External';

P0 = addoutput(niIO,settings.devID,'port0/line0','Digital');
P0.Name = 'Iontophoresis';
output(start_length*niIO.Rate + [0:stim_length*niIO.Rate],3) = 1;


%% record trial, save parameters          
trialData   = cell(trial_repeats, 1);                  %preallocate a cell to store the recording of each trial
rawData     = cell(trial_repeats, 1);

figure(1); clf;
h(1) = subplot(4,1,1:3);
h(2) = subplot(4,1,4);

figure(2); clf
h1 = subplot(3,1,1:2);     hold on;    ylabel('E (mV)');
                yyaxis(h1,'right');    ylabel('I (mV)');
h2 = subplot(3,1,3);       hold on;    ylabel('Input Resistance (M\Omega)');
                yyaxis(h2,'right');    ylabel('baseline potential (mV)');

t   = 0;
drugs = {};
while t < trial_repeats
    t = t+1;
    fprintf(['\n********** Acquiring Trial ', num2str(t),' ***********\n'])
    
    [rawDataTrial, trialTime]   = readwrite(niIO,output); %CHANGED FROM ELENAS CODE BECAUSE OF 2020 UPDATE
    rawData{t}                  = rawDataTrial; 
    
    [trialMeta.gain, trialMeta.mode, trialMeta.freq] = decodeTelegraphedOutput(rawDataTrial);
   
    % convert data into standard units with gain from ephysSettings
    trialData{t}.current = settings.current.softGain .* rawDataTrial.current; % pA
    trialData{t}.voltage = settings.voltage.softGain .* rawDataTrial.voltage; % mV
    trialData{t}.time    = rawDataTrial.Time;
    trialData{t}.output  = output;
    trialData{t}.datetime= datetime; 
    
    if inR_log
        [trialData{t}.InputData,trialData{t}.inR] = measureInputResistance(niIO,settings,2);
    else
        trialData{t}.inR = nan;
        trialData{t}.InputData = nan;
    end

    stim_start = find(diff(output(:,3)) > 0);
    stim_end   = find(diff(output(:,3)) < 0);
    fr         = 1e3;
    
    %% COPY AND PASTED PLOTTING STUFF
    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* rawDataTrial.s_output;  %mV
            
            % Plot vclamp trial
            plot(h(1),trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel(h(1),'Current (pA)')
            
            plot(h(2),trialData{t}.time, trialData{t}.voltage, 'k')
            ylabel(h(2),'Voltage (mV)')
            xlabel('Time (s)')
            
            figure(1)
            sgtitle(['Trial ' num2str(t)])
            linkaxes(h,'x')
                      
            smooth_data = medfilt1(trialData{t}.scaledOutput,fr);

            b = mean(smooth_data(1:stim_start)); %find mean of the data preceding the stim. window length of stim length
            v1 = max(smooth_data(stim_end:end)); %find mean of the data succeeding the stim. window length of stim length
            v2 = min(smooth_data(stim_end:end));
            
            yyaxis(h1,'left')
            scatter(h1,trialData{t}.datetime, (v1 - b),'filled')
            yyaxis(h1,'right')
            scatter(h1,trialData{t}.datetime, (v2 - b),'filled')
    
            yyaxis(h2,'left')
            scatter(h2,trialData{t}.datetime, max(trialData{t}.inR,0),'filled')
            yyaxis(h2,'right')
            scatter(h2,trialData{t}.datetime, mean(trialData{t}.scaledOutput),'filled')
            
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* rawDataTrial.s_output;  %mV
            
            
            % Plot iclamp trial
            plot(h(1),trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel(h(1),'Voltage (mV)')
            
            plot(h(2),trialData{t}.time, trialData{t}.current, 'k')
            ylabel(h(2),'Current (pA)')
            xlabel(h(2),'Time')
            
            figure(1)
            sgtitle(['Trial ' num2str(t)])
            linkaxes(h,'x')
            
            smooth_data = medfilt1(trialData{t}.scaledOutput,fr);
            b = mean(smooth_data(1:stim_start)); %find mean of the data preceding the stim. window length of stim length
            v1 = max(smooth_data(stim_end:end)); %find mean of the data succeeding the stim. window length of stim length
            v2 = min(smooth_data(stim_end:end));
            
            yyaxis(h1,'left')
            scatter(h1,trialData{t}.datetime, (v1 - b),'filled')
            yyaxis(h1,'right')
            scatter(h1,trialData{t}.datetime, (v2 - b),'filled')
            
            yyaxis(h2,'left')
            scatter(h2,trialData{t}.datetime, max(trialData{t}.inR,0),'filled')
            yyaxis(h2,'right')
            scatter(h2,trialData{t}.datetime, mean(trialData{t}.scaledOutput),'filled')
    end
    
    if t == trial_repeats
        if input('add drugs?: ')
            tmp = size(drugs,1) + 1;
            drugs{tmp,1} = input('name: ','s');
            drugs{tmp,2} = datetime;
        end
        
        tmp = input('add trials: ');
        trial_repeats = trial_repeats + tmp;
        trialData   = [trialData;cell(tmp, 1)];                  %preallocate a cell to store the recording of each trial
        rawData     = [rawData;cell(tmp, 1)];
    end
end
trialMeta.trialDuration_s = length(trialData{t}.time)/niIO.Rate;
trialMeta.trials      =  trial_repeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};
trialMeta.drugs       = drugs;

fprintf('\n********** acquireSimpleTrial Complete **********\n' )