function h = plot_traces(allData,output_channel)
h = figure;
subplot(3,1,1)
plot(allData.trialMeta.baseline{1}.time,allData.trialMeta.baseline{1}.scaledOutput,'k')
title('baseline')
set(gca,'FontSize',20)
    
idx = round(linspace(1,allData.trialMeta.trials,6));
for i = 1:6
    start_log       = diff(allData.trialData{idx(i)}.output(:,output_channel)) > 0;
    stim_start      = seconds(allData.trialData{idx(i)}.time(start_log));
    end_log         = diff(allData.trialData{idx(i)}.output(:,output_channel)) < 0;
    stim_end        = seconds(allData.trialData{idx(i)}.time(end_log));
    
    
    subplot(3,3,i+3)
    plot(allData.trialData{idx(i)}.time,allData.trialData{idx(i)}.scaledOutput,'k')
    y = ylim;
    patch([stim_start,stim_start,stim_end,stim_end],[-100,100,100,-100], [.1, 0, 1], 'edgecolor','none', 'facealpha',.1)
    ylim([y(1),y(2)])
    title(sprintf('%s',datestr(allData.trialData{idx(i)}.datetime, 'HH:MM:SS')))
    tmp  = xticks; xticks(tmp( mod(tmp,seconds(1)) == seconds(0)));
    set(gca,'FontSize',20)
end