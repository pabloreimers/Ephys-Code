function [inputData,inputResistance] = measureInputResistance(niIO, settings, num_out, External_Channel)
% finds the current recording mode (current or voltage clamp). then step
% either 5pA or 5mV, and do a little V=IR to find the resistance.

%%
fprintf('\n^^^^^^^^^ Acquiring Input Resistance ^^^^^^^^^\n' )

%Record for 1/5 of a second, no outputs, just to observe recording mode,
%and set appropriate external command gain
output = zeros(niIO.Rate/5,num_out);
[rawDataTrial, ~]   = readwrite(niIO,output);
[~,tmp_mode,~] = decodeTelegraphedOutput(rawDataTrial);

switch tmp_mode
    case {'I=0','I-Clamp Normal','I-Clamp Fast'}
        gain = settings.daq.current_extGain;
    case {'Track','V-Clamp'}
        gain = settings.daq.voltage_extGain;
end

% Set analog output command as a -5mV or pA step. The appropriate voltage
% sent to the axopatch is set by the gain, above. Everything should be in
% mV or pA, as defined in ephysSettings.
%The protocol here is 1s of nothing, and 1s of injection
output  = zeros(3*niIO.Rate,num_out);
output(1*niIO.Rate:2*niIO.Rate,External_Channel) = -5;
output  = output / gain; 

%% ACQUIRE TRIAL, CALCULATE INPUT RESISTANCE

[rawData, trialTime]   = readwrite(niIO,output);

[trialMeta.gain,trialMeta.mode,trialMeta.freq]= decodeTelegraphedOutput(rawData);

% Process non-scaled data, change raw data channels based on your setup
inputData(:,1) = settings.current.softGain .* rawData.current;
inputData(:,2) = settings.voltage.softGain .* rawData.voltage;
I = inputData(:,1);
V = inputData(:,2);

switch trialMeta.mode
    % Voltage Clamp
    case {'Track','V-Clamp'}
        settings.scaledOutput.softGain = 1e3 / (settings.current.beta * trialMeta.gain); %in V-clamp, I = alpha*beta mV / pA. So to get pA recorded, we turn the scaled output (in Volts) into mV and divide by alpha (trialMeta.gain, set in scaled output in axopatch) and beta (settings.current.beta, set in config on axopatch)
        inputData(:,3) = settings.scaledOutput.softGain .* rawData.s_output;  %mV
        I = inputData(:,3);
        % Current Clamp
    case {'I=0','I-Clamp Normal','I-Clamp Fast'}
        settings.scaledOutput.softGain = 1e3 / (trialMeta.gain); %in I-clamp, Vm = alpha mV per mV. alpha is the gain decoded in telegraphed output. to get mV recorded, we turn the scaled output (in Volts) into mV and divide by alpha (trialMeta.gain)
        inputData(:,3) = settings.scaledOutput.softGain .* rawData.s_output;  %pA
        V = inputData(:,3);
end

dI = mean(I(output(:,External_Channel) ~= 0)) - mean(I(output(:,External_Channel) == 0)); % mean high current value - mean low current value
dV = mean(V(output(:,External_Channel) ~= 0)) - mean(V(output(:,External_Channel) == 0)); % mean high voltage value - mean low voltage value

inputResistance = (dV * 1e-3) / (dI * 1e-12) * 1e-6;                           %convert mV to Volts, pA to Amps, and Ohms to megaOhms

fprintf(['\n^^^^^^^^^^^^^^ Rinput = ' ,num2str(round(inputResistance)),' MOhm ^^^^^^^^^^^^^\n'])
