function pipetteResistance = measurePipetteResistance()

fprintf('\n^^^^^^^^ Acquiring Pipette Resistance ^^^^^^^^\n' )

%take a sample measurement with seal test on, recording the current
%required to make 5mV steps (in votlage clamp)
[~, trialData, ~] = acquireSimpleTrial(1,.5,0,0,0);

%define current and voltage outside of structure for easy access
I = trialData{1}.scaledOutput;  %mV
V = trialData{1}.voltage;

%make a logical referencing when the voltage is at the high step
v_log   = V > mean(V);

dV = median(V(v_log)) - median(V(~v_log));
dI = median(I(v_log)) - median(I(~v_log));

%R = V/I, in standard units. x ohms = y volts/ z amps. current and voltage
%are measured in mV and pA, so apply appropriate multipliers to convert
%back to standard units. report resistance in megaohms.
pipetteResistance = (dV * 1e-3) / (dI * 1e-12) * 1e-6;