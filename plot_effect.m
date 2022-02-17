function h = plot_effect(allData,output_channel)
fr = 1e3;
h = figure;
start_log       = find(diff(allData.trialData{1}.output(:,output_channel)) > 0);
end_log         = find(diff(allData.trialData{1}.output(:,output_channel)) < 0);


h1 = subplot(3,1,1:2);     hold on;    ylabel('E (mV)');
                yyaxis(h1,'right');    ylabel('I (mV)');
h2 = subplot(3,1,3);       hold on;    ylabel('Input Resistance (M\Omega)');
                yyaxis(h2,'right');    ylabel('baseline potential (mV)');
window_date  = NaT(allData.trialMeta.trials,1);
window_data  = nan(allData.trialMeta.trials,1);

for i = 1:allData.trialMeta.trials
%     dI = mean(allData.trialData{i}.InputData(2e4:4e4,1)) - mean(allData.trialData{i}.InputData([1:2e4,4e4:6e4],1));
%     dV = mean(allData.trialData{i}.InputData(2e4:4e4,2)) - mean(allData.trialData{i}.InputData([1:2e4,4e4:6e4],2));
%     allData.trialData{i}.inR = (dV * 1e-3) / (dI * 1e-12) * 1e-6;
    
    smooth_data = medfilt1(allData.trialData{i}.scaledOutput,fr);
    
    b = mean(smooth_data(1:start_log)); %find mean of the data preceding the stim. window length of stim length
    v1 = max(smooth_data(end_log:end)); %find mean of the data succeeding the stim. window length of stim length
    v2 = min(smooth_data(end_log:end));
    
    den = 1; %abs(mean(allData.trialData{i}.scaledOutput));
    
    yyaxis(h1,'left')
    scatter(h1,allData.trialData{i}.datetime, (v1 - b)/den,'filled')
    yyaxis(h1,'right')
    scatter(h1,allData.trialData{i}.datetime, (v2 - b)/den,'filled')
    
    
    yyaxis(h2,'left')
    scatter(h2,allData.trialData{i}.datetime, max(allData.trialData{i}.inR,0),'filled')
    yyaxis(h2,'right')
    scatter(h2,allData.trialData{i}.datetime, mean(allData.trialData{i}.scaledOutput),'filled')
    
    window_date(i) = allData.trialData{i}.datetime;
    window_data(i) = v1-b;
    
end
yyaxis(h1,'left'); %plot(h1,window_date,medfilt1(window_data,5),'r-')
xlim([allData.trialData{1}.datetime,allData.trialData{end}.datetime])

%for j = 1:2
subplot(3,1,1:2)
yyaxis left
y = yticks;
y = linspace(y(1),y(end),20);
x = xlim;
for i = 1:size(allData.trialMeta.drugs,1)
    x(1) = allData.trialMeta.drugs{i,2};
    plot(x,[y(end-i+1),y(end-i+1)],'-k','linewidth',4)
    text(x(1),y(end-i+1),allData.trialMeta.drugs{i,1},'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',20)
end
%end
set(h1,'FontSize',20)
set(h2,'FontSize',20)
linkaxes([h1,h2],'x')
xlim([allData.trialData{1}.datetime,allData.trialData{end}.datetime])
shg