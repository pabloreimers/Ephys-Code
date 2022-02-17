function [rawData, trialData, trialMeta] = acquireRunningTrial(trial_repeats, start_length, stim_length, end_length, inR_log)
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
A0 = addoutput(niIO,settings.devID,'ao0','Voltage'); %add an output channel to our NI breakout board
A0.Name = 'Shutter';


%% record trial, save parameters

% create vector of analog outputs, the output will be written at each sample time for stim or no stim(however many seconds per however many samples per second)
output      = [zeros(ceil(start_length * niIO.Rate),1);...
                 5*ones(ceil(stim_length * niIO.Rate),1);
                 zeros(ceil(end_length * niIO.Rate),1)]; 

%if recording input resistance, add an external command output channel
if inR_log
    A1 = addoutput(niIO,settings.devID,'ao1','Voltage');
    A1.Name = 'External';
    output   = [output,zeros(size(output))];
end
           
trialData   = cell(trial_repeats, 1);                  %preallocate a cell to store the recording of each trial
rawData     = cell(trial_repeats, 1);

figure(1); clf;
h(1) = subplot(4,1,1:3);
h(2) = subplot(4,1,4);

figure(2); clf
hold on
h2(1) = gca;

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
        [trialData{t}.InputData,trialData{t}.inR] = measureInputResistance(niIO,settings,trialMeta.mode);
    else
        trialData{t}.inR = nan;
        trialData{t}.InputData = nan;
    end

    stim_start = find(diff(output(:,1)) > 0);
    stim_end   = find(diff(output(:,1)) < 0);
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
                      
            smooth_data = movmean(trialData{t}.scaledOutput,fr);
            b = smooth_data(stim_start);
            v = smooth_data(stim_end);%define effect size as window after stim that is length of stim
            
            yyaxis(h2(1),'left')
            scatter(h2(1),datetime, v-b , 'k','filled'); 
            xticks([])
            
            if inR_log
                yyaxis(h2(1),'right')
                scatter(h2(1),datetime, trialData{t}.inR ,'filled'); 
                ylabel(h2(1),'Input Resistance (M\Omega)')
            end
            
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
            
            figure(2);
            hold on
            smooth_data = movmean(trialData{t}.scaledOutput,fr);
            b = mean(smooth_data(2*stim_start-stim_end:stim_start));
            v = mean(smooth_data(stim_end:(2*stim_end-stim_start)));%define effect size as window after stim that is length of stim
            yyaxis(h2(1),'left')
            ylabel('effect size (mV)')
            scatter(h2(1),datetime, v-b , 'k','filled'); 
            xticks([])
            if inR_log
                yyaxis(h2(1),'right')
                scatter(h2(1),datetime, trialData{t}.inR ,'filled'); 
                ylabel('Input Resistance (M\Omega)')
            end
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