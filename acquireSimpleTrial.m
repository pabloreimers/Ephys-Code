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
A0.Name = 'Shutter';

A1 = addoutput(niIO,settings.devID,'ao1','Voltage'); %add an output channel to our NI breakout board
A1.Name = 'External command';

num_out = 2;            %how many total output channels do we have
ext_channel = 2;   %what is the index of the external command output channel (on the niDAQ)
stim_channel= 1;   %what is the index of the channel that we want to output during a trial

%% check input resistance
if inR_log
    [~, trialMeta.inputR] = measureInputResistance(niIO,settings,num_out,ext_channel);
else
  trialMeta.inputR = nan;
end

%% record trial, save parameters

% create vector of analog outputs, the output will be written at each sample time for stim or no stim(however many seconds per however many samples per second)
output      = zeros(niIO.Rate * (start_length+stim_length+end_length), num_out);
output(round(niIO.Rate*start_length) : round(niIO.Rate*(start_length+stim_length))-1, stim_channel) = 5; %after the start length, through the duration of the stim length, pass a 5V analog output
           
trialData   = cell(trial_repeats, 1);                  %preallocate a cell to store the recording of each trial
rawData     = cell(trial_repeats, 1);

for t = 1:trial_repeats
    fprintf(['\n********** Acquiring Trial ', num2str(t),' ***********\n'])
    
    [rawDataTrial, trialTime]   = readwrite(niIO,output);
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
            settings.scaledOutput.softGain = 1e3 / (settings.current.beta * trialMeta.gain); %in V-clamp, I = alpha*beta mV / pA. So to get pA recorded, we turn the scaled output (in Volts) into mV and divide by alpha (trialMeta.gain, set in scaled output in axopatch) and beta (settings.current.beta, set in config on axopatch)
            trialData{t}.scaledOutput =  settings.scaledOutput.softGain .* rawDataTrial.s_output;  %mV
            
            % Plot vclamp trial
            figure(1); clf;
            h(1) = subplot(4,1,1:3);
            plot(trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,4);
            plot(trialData{t}.time, trialData{t}.voltage, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            linkaxes(h,'x')
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1e3 / (trialMeta.gain); %in I-clamp, Vm = alpha mV per mV. alpha is the gain decoded in telegraphed output. to get mV recorded, we turn the scaled output (in Volts) into mV and divide by alpha (trialMeta.gain)
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