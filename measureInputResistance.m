function [inputData,inputResistance] = measureInputResistance(niIO, settings, External_Channel)
% Aquires a trial of current clamp data obtained, with the external command 
% on, and then use this response to calculate the input resistance. The 
% obtained trace and calculated input resistance are then saved

%%
fprintf('\n^^^^^^^^^ Acquiring Input Resistance ^^^^^^^^^\n' )

% Set analog output command
currentPulse = zeros(niIO.Rate,1);
%tmp = [1:2000:niIO.Rate] + [1:1000]';
%tmp = [1:ceil(niIO.Rate/2):niIO.Rate] + [1:ceil(niIO.Rate/4)]'; %do a pulse at 4 hertz
%tmp = tmp(:);
%currentPulse(tmp) = 1;
currentPulse = [ones(niIO.Rate,1);zeros(niIO.Rate,1)];                          % reach steady state
currentPulse = currentPulse * -.5;                                              % convert stim protocol to pA
currentPulse = currentPulse * (1/settings.axopatch_picoAmps_per_volt);          % convert pA to Vout
%currentPulse = currentPulse * (1/settings.AO_output_scaling_factor); % not
%sure why elena had this line in here, but it messes mine up
% if ~strcmp(mode,'I-Clamp Normal')
%     currentPulse = currentPulse * 100;
% end
output = zeros(size(currentPulse,1),3); % initialize
output(:,External_Channel) = currentPulse;

%% ACQUIRE TRIAL, CALCULATE INPUT RESISTANCE

[rawData, trialTime]   = readwrite(niIO,output);

[trialMeta.gain,trialMeta.mode,trialMeta.freq]= decodeTelegraphedOutput(rawData);

% Process non-scaled data, change raw data channels based on your setup
inputData(:,1) = settings.current.softGain .* rawData.current;
inputData(:,2) = settings.voltage.softGain .* rawData.voltage;
I = settings.current.softGain .* rawData.current;
V = settings.voltage.softGain .* rawData.voltage;

switch trialMeta.mode
    % Voltage Clamp
    case {'Track','V-Clamp'}
        settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
        inputData(:,3) = settings.scaledOutput.softGain .* rawData.s_output;  %mV
        I = inputData(:,3);
        % Current Clamp
    case {'I=0','I-Clamp Normal','I-Clamp Fast'}
        settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
        inputData(:,3) = settings.scaledOutput.softGain .* rawData.s_output;  %pA
        V = inputData(:,3);
end

% iDelta = mean(inputData((1*niIO.Rate):(2*niIO.Rate),1)) - mean(inputData(1:(2*niIO.Rate+1),1));
% vDelta = mean(inputData((1*niIO.Rate):(6*niIO.Rate),3)) - mean(inputData(1:(2*niIO.Rate+1),3));
dI = mean(I(output(:,External_Channel) ~= 0)) - mean(I(output(:,External_Channel) == 0)); % mean high current value - mean low current value
dV = mean(V(output(:,External_Channel) ~= 0)) - mean(V(output(:,External_Channel) == 0)); % mean high voltage value - mean low voltage value

inputResistance = (dV * 1e-3) / (dI * 1e-12) * 1e-6;                           %convert mV to Volts, pA to Amps, and Ohms to megaOhms

fprintf(['\n^^^^^^^^^^^^^^ Rinput = ' ,num2str(round(inputResistance)),' MOhm ^^^^^^^^^^^^^\n'])
