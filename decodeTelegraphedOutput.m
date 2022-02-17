function [gain,mode,freq] = decodeTelegraphedOutput(data,gainCh,modeCh,freqCh)
% decodeTelegraphedOutput decodes amplifier input contianing the mode,
% gain, and frequency
ephysSettings;
gainCh = settings.raw.gain;
modeCh = settings.raw.mode;
freqCh = settings.raw.frequency;


%% Gain
gainOut = mean(data.gain);

if gainOut> 0        && gainOut< 0.75
    gain = 0.05;
elseif gainOut> 0.75 && gainOut< .80
    gain = 0.1;
elseif gainOut> .80 && gainOut< 1.00
    gain = 0.2;
elseif gainOut> 1.50 && gainOut< 2.00
    gain = 0.5;
elseif gainOut> 2.00 && gainOut< 2.50
    gain = 1;
elseif gainOut> 2.50 && gainOut< 3.00
    gain = 2;
elseif gainOut> 3.00 && gainOut< 3.60
    gain = 5;
elseif gainOut> 3.60 && gainOut< 4.20
    gain = 10;
elseif gainOut> 4.20 && gainOut< 4.70
    gain = 20;
elseif gainOut> 4.70 && gainOut< 5.20
    gain = 50;
elseif gainOut> 5.20 && gainOut< 5.80
    gain = 100;
elseif gainOut> 5.80 && gainOut< 6.50
    gain = 200;
elseif gainOut> 6.50 && gainOut< 7.00
    gain = 500;
else
    gain = Inf;
end

%% freq 
freqOut  = mean(data.frequency);

if freqOut > 0 && freqOut < 3
    freq = 1;
elseif freqOut > 3 && freqOut < 5
    freq = 2;
elseif freqOut > 5 && freqOut < 7
    freq = 5;
elseif freqOut > 7 && freqOut < 9
    freq = 10;
elseif freqOut > 9
    freq = 100;
end

%% Mode
modeOut = mean(data.mode);

if modeOut> 0 && modeOut< 1.5
    mode = 'I-Clamp Fast';
elseif modeOut> 1.5 && modeOut< 2.5
    mode = 'I-Clamp Normal';
elseif modeOut> 2.5 && modeOut< 3.5
    mode = 'I=0';
elseif modeOut> 3.5 && modeOut< 5
    mode = 'Track';
elseif modeOut> 5
    mode = 'V-Clamp';
end
