%% ephysSettings
%specifies hard-coded parameters for ephys acquisition includig data
%directory path and axopatch parameters
% PR 2021/21/1

%% Device Configuration
settings    = struct();         %initialize a struct to store all settings

comptype    = computer;         %store the kind of computer

settings.mainDataDir = 'C:\Users\ReimersPabloAlejandr\Documents\Code\pablo\Data'; %set path to store data

settings.devID      = 'Dev1';   % NO IDEA WHAT THIS DOES
settings.sampRate   = 2e4;

%% NI breakout board channel assignments
settings.bob.current    = 0;
settings.bob.voltage    = 1;
settings.bob.s_output   = 2;
settings.bob.mode       = 3;
settings.bob.gain       = 4;
settings.bob.frequency  = 5;

settings.bob.channelIDs = 0:5;  %FLAG THIS FOR ELENA, MIGHT BE OUT OF ORDER

%% Format of the raw data. used in DAQsettings.m
settings.raw.current    = 1;
settings.raw.voltage    = 2;
settings.raw.s_output   = 3;
settings.raw.mode       = 4;
settings.raw.gain       = 5;
settings.raw.frequency  = 6;

%% Current and voltage settings
settings.current.beta       = 1; %this is defined in config. on the front of the axopatch
settings.current.gain       = 1 * settings.current.beta; %rear switch for current output (I 1kHz) set to beta * mV/pA. 

settings.current.softGain   = 1e3*(1/(settings.current.gain)); % mV/V * (pA/mV) = pA / V
%this number is a multiplier on the voltage coming out of the axopatch amplifier in the
%current channel (in the back), so that it can be appropriately converted to the pA
%recorded by the axopatch amplifier

settings.voltage.gain       = 10; % 10 Vm, set coming out of the back of the amplifier. This is just labelled, no way to change.
settings.voltage.softGain   = 1e3 / ( settings.voltage.gain); % To get voltage in mV
%this number is a multiplier on teh voltage coming out of he axopatch
%amplifier in the 10Vm Voltage channel (in the back), so that it can be
%appropriately converted to the mV recorded by the axopatch amplifier.

%COPY AND PASTE BELOW. MAKE SURE TO GET CLARITY
% Digital Voltage output settings:
settings.daq.current_extGain = 2 / settings.current.beta * 1e3; % external command of injecting (2/beta) nA/V sent to axopatch. This is the same number as pA / mV sent to axopathc.  
%if I want to inject 5pA, I need to send (5 pA / current_extGain pA/mV) to the
%axopatch. I AM SENDING MV TO AXOPATCH.
settings.daq.voltage_extGain = 20e-3 / 1 * 1e3; %external commend of inject (20mV/V) sent to axopatch. For every V sent to axopatch, the holding current jumps by 20mV. 
%if I want to step the holding voltage by 5mV, I need to send 1/4 of a V, 250mV, or 5mV / voltage_extGain
%to the axopatch. I AM SENDING MV TO AXOPATCH


%% OLD CODE TRYING TO OVERWRITE
% settings.daq.frontExtScale = 20 / 1000; %20mV/ 1000mV (1V) amplifier cuts the voltage down by this factor, every 1volt from the DAQ is 20mV into the Axopatch
% settings.daq.voltageConversionFactor =  1 / (settings.daq.frontExtScale * settings.daq.voltageDividerScaling * 1e3); % use this for votlage clamp experiment commands 
% %1 Volt = 2nA * Beta (1 normally)
% 
% % voltage clamp mode -
% %   "20 mV/V", +1 V input produces +20 mV
% %   input voltage range of -10 to +10 V produces -1000 to +1000 mV
% settings.axopatch_mV_per_volt=20/1; % 20 mV / 1 V 
% 
% % current clamp mode -
% %    "2 / (beta nA/V)" & beta = 1 nA/V
% %    +1 V produces 2 nA (2000 pA)
% %    input voltage range of -10 to +10 V produces -20 to +20 nA (-20,000 to +20,000 pA) 
% settings.axopatch_picoAmps_per_volt=2000/1;
% 
% 
% % If want to use -10 to +10 V (AO command) to inject -100 to +100 pA of current, 
% % we need to scale the output of the daq board by a factor of 100/20,000=0.0050.
% % voltage divider configuration is here:
% % http://www.falstad.com/circuit/circuitjs.html?cct=$+1+0.000005+10.20027730826997+63+10+62%0Ar+288+336+288+256+0+4700%0Aw+288+256+368+160+0%0Aw+288+256+288+160+0%0Aw+288+96+288+64+0%0Aw+160+336+160+272+0%0Ag+432+96+464+64+0%0Aw+288+160+288+96+0%0AR+160+272+160+224+0+0+40+1+0+0+0.5%0Aw+160+336+160+400+0%0Aw+160+400+288+400+0%0Ar+288+400+288+336+0+4700%0Ar+368+160+400+128+0+22%0Ar+400+128+432+96+0+22%0A
% settings.AO_output_scaling_factor=0.00476;
