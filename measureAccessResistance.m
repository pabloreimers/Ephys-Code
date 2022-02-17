function accessResistance = measureAccessResistance()

fprintf('\n^^^^^^^^ Acquiring Access Resistance ^^^^^^^^\n' )

%take a sample measurement with seal test on, recording the current
%required to make 5mV steps (in votlage clamp)
[~, trialData, ~] = acquireSimpleTrial(1,.5,0,0,0);

%define current and voltage outside of structure for easy access
I = trialData{1}.current;
V = trialData{1}.voltage;

%make a logical referencing when the voltage is at the high step
v_log   = V > mean(V);
up_dist  = min(diff(find(diff(v_log) > 0))); %make logical of start and end of up and down. extract just the up indexes. find min distance between ups. 

peakCurrent = median(findpeaks(I,'minpeakdistance',up_dist)); %find the median peak current. 
baseCurrent = median(I(~v_log)); %find the current before the voltage pulse. median to be robust to 1) noise and 2) capacitance


dV = median(V(v_log)) - median(V(~v_log));
dI = peakCurrent - baseCurrent;

%R = V/I, in standard units. x ohms = y volts/ z amps. current and voltage
%are measured in mV and pA, so apply appropriate multipliers to convert
%back to standard units. report resistance in megaohms.
accessResistance = (dV * 1e-3) / (dI * 1e-12) * 1e-6;