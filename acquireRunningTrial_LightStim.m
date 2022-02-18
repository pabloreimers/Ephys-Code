function [rawData, trialData, trialMeta] = acquireRunningTrial_LightStim(trial_repeats, start_length, stim_length, end_length)
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
[~, trialMeta.inputR_start] = measureInputResistance(niIO,settings,num_out,ext_channel);

%% Initialize analog output matrix, experimental parameters, and figures

% create vector of analog outputs, the output will be written at each sample time for stim or no stim(however many seconds per however many samples per second)
output      = zeros(niIO.Rate * (start_length+stim_length+end_length), num_out);
stim_idx    = round(niIO.Rate*start_length : niIO.Rate*(start_length+stim_length));
output(stim_idx, stim_channel) = 5; %after the start length, through the duration of the stim length, pass a 5V analog output
           
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

%% Acquire trials
while t < trial_repeats
    t = t+1;
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
    
    [trialData{t}.InputData,trialData{t}.inR] = measureInputResistance(niIO,settings,num_out,ext_channel);

    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1e3 / (settings.current.beta * trialMeta.gain); %in V-clamp, I = alpha*beta mV / pA. So to get pA recorded, we turn the scaled output (in Volts) into mV and divide by alpha (trialMeta.gain, set in scaled output in axopatch) and beta (settings.current.beta, set in config on axopatch)
            trialData{t}.scaledOutput =  settings.scaledOutput.softGain .* rawDataTrial.s_output;  %mV
            
            % Plot vclamp trial data
            plot(h(1),trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel(h(1),'Current (pA)')
            
            plot(h(2),trialData{t}.time, trialData{t}.voltage, 'k')
            ylabel(h(2),'Voltage (mV)')
            xlabel('Time (s)')
            
            figure(1)
            sgtitle(['Trial ' num2str(t)])
            linkaxes(h,'x')          
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1e3 / (trialMeta.gain); %in I-clamp, Vm = alpha mV per mV. alpha is the gain decoded in telegraphed output. to get mV recorded, we turn the scaled output (in Volts) into mV and divide by alpha (trialMeta.gain)
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
    end
    
            %plot online analysis
            
            %cell stats (input resistance and baseline potential)
            yyaxis(h2,'left')
            scatter(h2,trialData{t}.datetime, trialData{t}.inR,'filled')
            yyaxis(h2,'right')
            scatter(h2,trialData{t}.datetime, mean(trialData{t}.voltage),'filled')
                
            %plot min and max of median filtered (smoothed) traces
            %(baseline subtracted)
            stim_start = find(diff(output(:,stim_channel)) > 0);
            stim_end   = find(diff(output(:,stim_channel)) < 0);
            fr         = 1e3;
            smooth_data = medfilt1(trialData{t}.scaledOutput,fr);
            b  = mean(smooth_data(1:stim_start)); %find mean of the data preceding the stim. window length of stim length
            v1 = max(smooth_data(stim_end:end)); %find mean of the data succeeding the stim. window length of stim length
            v2 = min(smooth_data(stim_end:end));
            
            yyaxis(h1,'left')
            scatter(h1,trialData{t}.datetime, (v1 - b),'filled')
            yyaxis(h1,'right')
            scatter(h1,trialData{t}.datetime, (v2 - b),'filled')
    
   %user input add trials
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

fprintf('\n********** acquireRunningTrial_LightStim Complete **********\n' )