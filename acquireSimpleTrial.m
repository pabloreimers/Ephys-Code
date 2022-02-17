function [rawData, trialData, trialMeta] = acquireSimpleTrial(trial_repeats, start_length, stim_length, end_length, inR_log)
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
A0.Name = 'External command';

%% check input resistance
if inR_log
    [~, trialMeta.inputR] = measureInputResistance(niIO,settings,1);
else
  trialMeta.inputR = nan;
end

%% record trial, save parameters

% create vector of analog outputs, the output will be written at each sample time for stim or no stim(however many seconds per however many samples per second)
output      = [zeros(ceil(start_length * niIO.Rate),1);...
               5*ones(ceil(stim_length * niIO.Rate),1);
               zeros(ceil(end_length * niIO.Rate),1)]; 
           
trialData   = cell(trial_repeats, 1);                  %preallocate a cell to store the recording of each trial
rawData     = cell(trial_repeats, 1);

for t = 1:trial_repeats
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
    
    
    %% COPY AND PASTED PLOTTING STUFF
    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* rawDataTrial.s_output;  %mV
            
            % Plot vclamp trial
            figure(1); clf;
            h(1) = subplot(4,1,1:3);
            plot(trialData{t}.time, trialData{t}.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,4);
            plot(trialData{t}.time, trialData{t}.voltage, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            linkaxes(h,'x')
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* rawDataTrial.s_output;  %mV
            
            % Plot iclamp trial
            figure(1); clf;
            h(1) = subplot(4,1,1);
            plot(trialData{t}.time, trialData{t}.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,2:4);
            plot(trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            sgtitle(['Trial ' num2str(t)])
            
            linkaxes(h,'x')
    end
end
trialMeta.trialDuration_s = length(trialData{t}.time)/niIO.Rate;
trialMeta.trials      =  trial_repeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};

fprintf('\n********** acquireSimpleTrial Complete **********\n' )