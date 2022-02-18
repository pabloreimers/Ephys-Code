%% Generate DAQ Sessions and Establish Input/Output channels

fprintf('\n********** Setting Parameters ***********\n')

niIO    = daq('ni');    %initialize input/output object for data aqcuisition board. from matlab data acquisition toolbox
devID   = settings.devID;
niIO.Rate = settings.sampRate;      %sample rate in Hz

%prepare analog input channels according to DAC break out board
AI              = addinput(niIO,devID,settings.bob.channelIDs,'Voltage'); %add inputs to our input output object. they will expect analog input in the form of a voltage
channel_names   = fieldnames(settings.raw); %extract the channel names from the settings struct

for i = 1:length(channel_names) %for each channel, store its name.
    AI(i).Name = channel_names{i};
    AI(i).TerminalConfig = 'SingleEnded'; %no idea what this does
end
