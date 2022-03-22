%% collect experiment titles (likely drugs or experimental differences) and preprocess the data

k           = input('Number of experiments:  ');
exp_cell   = cell(2,k);

for i = 1:k
    %exp_cell{1,i}   = input('Experiment Name:   ','s');
    fprintf('select experiment data\n')
    [f,p]           = uigetfile('.mat','Select Data', 'C:\Users\ReimersPabloAlejandr\Documents\Code\pablo\Data\');
    tmp             = regexp(p,'\','split');
    exp_cell{1,i}   = tmp{end-2};
    exp_cell{2,i}   = load([p,f]);
    n               = length(exp_cell{2,i}.allData.trialData);
    tmp             = nan(n, length(exp_cell{2,i}.allData.trialData{1}.voltage),5);
    for ii = 1:n
        tmp(ii,:,1) = exp_cell{2,i}.allData.trialData{ii}.voltage;
        tmp(ii,:,2) = exp_cell{2,i}.allData.trialData{ii}.current;
        tmp(ii,:,3) = seconds(exp_cell{2,i}.allData.trialData{ii}.time);
        tmp(ii,:,4) = exp_cell{2,i}.allData.trialData{ii}.output(:,find(any(exp_cell{2,i}.allData.trialData{ii}.output,1)));
        tmp(ii,:,5) = exp_cell{2,i}.allData.trialData{ii}.scaledOutput;
    end
    D.voltage       = tmp(:,:,1);
    D.current       = tmp(:,:,2);
    D.time          = tmp(:,:,3);
    D.output        = tmp(:,:,4);
    D.scaledOutput  = tmp(:,:,5);
    
    exp_cell{3,i}   = D;
end



%% plot all trials, mean and s.e.m show traces and heatmap baseline subtracted
max_effect          = nan(k,2);
fr                  = allData.trialMeta.daqRate * (5e-3);              %smoothing window size. data captured at 20e4 frames/s, so 20e4 f / 1000 ms, so (20e4/1000) * 5 frames in a 5ms window

figure(1); clf
figure(2); clf
for i = 1:k
    t               = exp_cell{3,i}.time(1,:);
    start_log       = diff(any(exp_cell{3,i}.output,1)) > 0;
    stim_start      = t(start_log);
    end_log         = diff(any(exp_cell{3,i}.output,1)) < 0;
    stim_end        = t(end_log);
    stim_frames     = sum(exp_cell{3,i}.output(1,:)>0);
    v               = medfilt1(exp_cell{3,i}.voltage,fr,[],2,'truncate');          %smooth the data with a moving window of size fr. 
    sem             = std(v,1) / sqrt(size(v,1));
    b               = mean(v(:, 1:find(start_log)),2);
    n               = size(v,1);
    mean_curr       = mean(reshape(exp_cell{3,i}.current,1,[]));
    
    figure(1)
    subplot(3,k,[i,i+k]);
    hold on
    title(exp_cell{1,i})
    plot(t, v - b, 'color', [0,0,0, 1/size(exp_cell{3,i}.voltage,1)])
    plot(t, mean(v - b,1), 'r')
    patch([t,fliplr(t)], [mean(v - b,1) + sem, fliplr(mean(v - b,1) - sem)],[1,0,0],'edgecolor','none','facealpha',0.5)
    y = ylim;
    patch([stim_start,stim_start,stim_end,stim_end],[y(1),y(2),y(2),y(1)], [.1, 0, 1], 'edgecolor','none', 'facealpha',.1)
    ylim(y)
    x = xlim;
    xlabel('time (s)')
    text(x(1) + .1*(x(2)-x(1)),y(1) + .1*(y(2)-y(1)),sprintf('rest: %.1f mV\nI: %.1f pA',mean(b),mean_curr))
    axis tight
    subplot(3,k,i+2*k)
    imagesc(-(max(v- b,0)))
    colormap('bone')
    xticks([])
    
    figure(2)
    hold on
    smooth_v = v;
    max_idx = find(end_log)+[1:stim_frames];
    
    max_effect(i,1) = mean(max(smooth_v(:,max_idx),[],2) - b); %define the max effect size as the mean subtracted membane voltage at the offset of the stim
    max_effect(i,2) = std(max(smooth_v(:,max_idx),[],2) - b)/sqrt(size(v,1));
    scatter(i*ones(n,1) + [1/n:1/n:1]' -0.5,  max(smooth_v(:,max_idx),[],2) - b,20,'MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',0.25,'MarkerEdgeColor','none')
    errorbar(i, max_effect(i,1),max_effect(i,2),'ko','MarkerSize',5,'MarkerFaceColor','red','MarkerEdgeColor','none')
end
figure(1)
subplot(3,k,[1,k+1])
ylabel('V_m (mV, baseSub)')
subplot(3,k,2*k+1)
ylabel('trial')

figure(2)
ylabel('effect size (mV)')
xticks(1:k)
xticklabels(exp_cell(1,:))

%% Plot by Current Injection
[f,p]           = uigetfile('.mat','Select Current Injection Data', 'C:\Users\ReimersPabloAlejandr\Documents\Code\pablo\Data\');
load([p,f]);

idx = any(allData.trialData{1}.output,1);
D.voltage = cell2mat(cellfun(@(s)(s.voltage),allData.trialData','UniformOutput',false))';
D.current = cell2mat(cellfun(@(s)(s.current),allData.trialData','UniformOutput',false))';
D.scaledOutput = cell2mat(cellfun(@(s)(s.scaledOutput),allData.trialData','UniformOutput',false))';
D.output  = cell2mat(cellfun(@(s)(s.output(:,idx)),allData.trialData','UniformOutput',false))';
%D.I_inj   = cell2mat(cellfun(@(s)(s.I_inj),allData.trialData','UniformOutput',false))';
D.time    = allData.trialData{1}.time';
D.inR     = cell2mat(cellfun(@(s)(s.inR),allData.trialData','UniformOutput',false))';

g = round(mean(D.current(:,logical(D.output(1,:))),2) - mean(D.current(:,~logical(D.output(1,:))),2));
%g = I_inj;
cmap = turbo(max(g) - min(g) + 1);
figure(1); clf; hold on
for i = unique(g)'
    idx = ceil(i - min(g))+1;
    plot(D.time,nan(size(D.time)),'Color',cmap(idx,:),'LineWidth',2)
end
legend([num2str(unique(g)),repmat(' pA',size(unique(g)))],'AutoUpdate','off')

for i = 1:length(g)
    idx = ceil(g(i) - min(g))+1;
    plot(D.time,D.voltage(i,:),'Color',cmap(idx,:))
end
ylabel('Membrane Voltage (mV)')
xlabel('time')


%% Plot individual traces
n_trace = 20;
x = ceil(sqrt(n_trace));
y = floor(sqrt(n_trace));

for i = 1:size(exp_cell,2)
    figure   
    n_trace = size(exp_cell{3,i}.voltage,1);
    for ii = 1:n_trace
        
        subplot(ceil(1.5*y), x, 2*x + ii)
        
        yyaxis right
        plot(exp_cell{3,i}.voltage(ii,:),'Color',[1,0.5,0,0.5])
        set(gca,'YColor',[1,0.5,0,0.5])
        
        yyaxis left
        plot(exp_cell{3,i}.current(ii,:),'Color',[0,0.5,1,1])
        set(gca,'YColor',[0,0.5,1,1])
        xticks([])
        axis tight
        

        %yticks([])
        
        if mod(ii, floor(n_trace/3)) == 0
            subplot(3,3,round(n_trace/ii))
            hold on
            plot(exp_cell{3,i}.voltage(ii,:),'Color',[.5,0.5,.5,0.5])
            plot(movmean(exp_cell{3,i}.voltage(ii,:),fr),'Color',[1,0.5,0])
            xticks([])
            ylabel('Vm (mV)')
            axis tight
        end
    end
        set(gcf,'Name',exp_cell{1,i})
end

%% plot baselines
k           = input('Number of experiments:  ');
exp_cell   = cell(2,k);
figure;

for i = 1:k
    %exp_cell{1,i}   = input('Experiment Name:   ','s');
    h3(i) = subplot(k,1,i);
    fprintf('select experiment data\n')
    [f,p]           = uigetfile();
    load([p,f])
    tmp             = regexp(p,'\','split');
    plot(allData.trialData{1}.time,allData.trialData{1}.voltage)
    hold on
    plot(allData.trialData{1}.time,movmean(allData.trialData{1}.voltage,fr))
    ylabel(tmp{end-2})
end

linkaxes(h3)
shg

%% load in running trial data
[f,p]           = uigetfile();
load([p,f]);

%% Plot eveerything (c-clamp)
figure;
subplot(3,1,1)
plot(allData.trialMeta.baseline{1}.time,allData.trialMeta.baseline{1}.scaledOutput,'k')
title('baseline')
set(gca,'FontSize',20)
    
idx = round(linspace(1,allData.trialMeta.trials,6));
for i = 1:6
    start_log       = diff(allData.trialData{idx(i)}.output(:,1)) > 0;
    stim_start      = seconds(allData.trialData{idx(i)}.time(start_log));
    end_log         = diff(allData.trialData{idx(i)}.output(:,1)) < 0;
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


figure;
h1 = subplot(3,1,1:2);     hold on;    ylabel('effect size (mV)');
                yyaxis(h1,'right');    ylabel('resting potential (mV)');
h2 = subplot(3,1,3);       hold on;    ylabel('Input Resistance (M\Omega)');
smooth_date  = NaT(allData.trialMeta.trials,1);
smooth_data  = nan(allData.trialMeta.trials,1);

for i = 1:allData.trialMeta.trials
    stim_start = find(start_log);
    stim_end   = find(end_log);
    b = mean(allData.trialData{i}.scaledOutput(2*stim_start - stim_end : stim_start)); %find mean of the data preceding the stim. window length of stim length
    v = mean(allData.trialData{i}.scaledOutput(stim_end : 2*stim_end - stim_start)); %find mean of the data succeeding the stim. window length of stim length
    
    yyaxis(h1,'left')
    scatter(h1,allData.trialData{i}.datetime, v - b,'filled')
    yyaxis(h1,'right')
    scatter(h1,allData.trialData{i}.datetime, mean(allData.trialData{i}.scaledOutput),'filled')
    
    scatter(h2,allData.trialData{i}.datetime, max(allData.trialData{i}.inR,0),'filled','k')
    
    smooth_date(i) = allData.trialData{i}.datetime;
    smooth_data(i) = v-b;
    
end
yyaxis(h1,'left'); plot(h1,smooth_date,movmean(smooth_data,5),'r-')

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
shg

%% plot current clamp 2
fr = 1e3;
figure;
subplot(3,1,1)
plot(allData.trialMeta.baseline{1}.time,allData.trialMeta.baseline{1}.scaledOutput,'k')
title('baseline')
ylabel('voltage (mV)')
set(gca,'FontSize',20)
    
idx = round(linspace(1,allData.trialMeta.trials,6));
for i = 1:6
    dim             = any(allData.trialData{idx(i)}.output,1);
    start_log       = diff(allData.trialData{idx(i)}.output(:,dim)) > 0;
    stim_start      = seconds(allData.trialData{idx(i)}.time(start_log));
    end_log         = diff(allData.trialData{idx(i)}.output(:,dim)) < 0;
    stim_end        = seconds(allData.trialData{idx(i)}.time(end_log));
    
    
    subplot(3,3,i+3)
    plot(allData.trialData{idx(i)}.time,allData.trialData{idx(i)}.scaledOutput,'k')
    y = ylim;
    patch([stim_start,stim_start,stim_end,stim_end],[y(1),y(2),y(2),y(1)], [.1, 0, 1], 'edgecolor','none', 'facealpha',.1)
    ylim([y(1),y(2)])
    title(sprintf('%s',datestr(allData.trialData{idx(i)}.datetime, 'HH:MM:SS')))
    tmp  = xticks; xticks(tmp( mod(tmp,seconds(1)) == seconds(0)));
    set(gca,'FontSize',20)
end


figure;
h1 = subplot(3,1,1:2);     hold on;    ylabel('effect size (mV)');
                yyaxis(h1,'right');    ylabel('baseline potential (mV)');
h2 = subplot(3,1,3);       hold on;    ylabel('Input Resistance (M\Omega)');
window_date  = NaT(allData.trialMeta.trials,1);
window_data  = nan(allData.trialMeta.trials,1);

for i = 1:allData.trialMeta.trials
%     dI = mean(allData.trialData{i}.InputData(2e4:4e4,1)) - mean(allData.trialData{i}.InputData([1:2e4,4e4:6e4],1));
%     dV = mean(allData.trialData{i}.InputData(2e4:4e4,2)) - mean(allData.trialData{i}.InputData([1:2e4,4e4:6e4],2));
%     allData.trialData{i}.inR = (dV * 1e-3) / (dI * 1e-12) * 1e-6;
    
    smooth_data = medfilt1(allData.trialData{i}.scaledOutput,fr);
    
    b = smooth_data(start_log); %find mean of the data preceding the stim. window length of stim length
    v = smooth_data(end_log); %find mean of the data succeeding the stim. window length of stim length
    den = 1; %abs(mean(allData.trialData{i}.scaledOutput));
    
    yyaxis(h1,'left')
    scatter(h1,allData.trialData{i}.datetime, (v - b)/den,'filled')
    yyaxis(h1,'right')
    scatter(h1,allData.trialData{i}.datetime, mean(allData.trialData{i}.scaledOutput),'filled')
    
    scatter(h2,allData.trialData{i}.datetime, max(allData.trialData{i}.inR,0),'filled','k')
    
    window_date(i) = allData.trialData{i}.datetime;
    window_data(i) = v-b;
    
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


%% Plot voltage clamp
fr = 1e3;
figure;
subplot(3,1,1)
%plot(allData.trialMeta.baselineI0.trialData{1}.time,allData.trialMeta.baselineI0.trialData{1}.scaledOutput,'k')
title('baseline')
ylabel('voltage (Vm)')
set(gca,'FontSize',20)
    
idx = round(linspace(1,allData.trialMeta.trials,6));
for i = 1:6
    start_log       = diff(allData.trialData{idx(i)}.output(:,1)) > 0;
    stim_start      = seconds(allData.trialData{idx(i)}.time(start_log));
    end_log         = diff(allData.trialData{idx(i)}.output(:,1)) < 0;
    stim_end        = seconds(allData.trialData{idx(i)}.time(end_log));
    
    
    subplot(3,3,i+3)
    plot(allData.trialData{idx(i)}.time,allData.trialData{idx(i)}.scaledOutput,'k')
    y = ylim;
    patch([stim_start,stim_start,stim_end,stim_end],[y(1),y(2),y(2),y(1)], [.1, 0, 1], 'edgecolor','none', 'facealpha',.1)
    ylim([y(1),y(2)])
    title(sprintf('%s',datestr(allData.trialData{idx(i)}.datetime, 'HH:MM:SS')))
    tmp  = xticks; xticks(tmp( mod(tmp,seconds(1)) == seconds(0)));
    set(gca,'FontSize',20)
end
subplot(3,3,4)
ylabel('current (pA)')


figure;
h1 = subplot(3,1,1:2);     hold on;    ylabel('effect size (pA)');
                yyaxis(h1,'right');    ylabel('baseline current (pA)');
h2 = subplot(3,1,3);       hold on;    ylabel('Input Resistance (M\Omega)');
window_date  = NaT(allData.trialMeta.trials,1);
window_data  = nan(allData.trialMeta.trials,1);
tmp = nan(allData.trialMeta.trials,3);
tmp2 = NaT(allData.trialMeta.trials,1);

for i = 1:allData.trialMeta.trials
    dI = median(allData.trialData{i}.InputData(20000:40000,3)) - median(allData.trialData{i}.InputData([1:20000,40000:60000],3)) ./ median(allData.trialData{i}.InputData(1:2e4,3));
    dV = median(allData.trialData{i}.InputData(20000:40000,2)) - median(allData.trialData{i}.InputData([1:20000,40000:60000],2));
    allData.trialData{i}.inR_scaled = (dV * 1e-3) / (dI * 1e-12) * 1e-6;
    
    smooth_data = medfilt1(allData.trialData{i}.scaledOutput,fr);
    
    b = smooth_data(start_log); %find mean of the data preceding the stim. window length of stim length
    v = smooth_data(end_log); %find mean of the data succeeding the stim. window length of stim length
    den = 1; %abs(mean(allData.trialData{i}.scaledOutput));
    
    tmp(i,1) = (v-b)/den;
    tmp(i,2) = mean(allData.trialData{i}.scaledOutput);
    tmp(i,3) = max(allData.trialData{i}.inR,0);
    tmp2(i) = allData.trialData{i}.datetime;
    
%     yyaxis(h1,'left')
%     scatter(h1,allData.trialData{i}.datetime, (v - b)/den,'filled')
%     yyaxis(h1,'right')
%     scatter(h1,allData.trialData{i}.datetime, mean(allData.trialData{i}.scaledOutput),'filled')
%     
%     scatter(h2,allData.trialData{i}.datetime, max(allData.trialData{i}.inR,0),'filled','k')
%     
    window_date(i) = allData.trialData{i}.datetime;
    window_data(i) = v-b;
    
end
yyaxis(h1,'left'); %plot(h1,window_date,medfilt1(window_data,5),'r-')
tmp(isoutlier(tmp,1)) = nan;

yyaxis(h1,'left'); scatter(h1,tmp2,tmp(:,1),'filled'); ylabel('effect size (pA)')
yyaxis(h1, 'right'); scatter(h1, tmp2,tmp(:,2),'filled'); ylabel('baseline current (pA');
scatter(h2,tmp2,tmp(:,3),'filled','k'); ylabel(h2,'Input Resistance (M\Omega)');
x = [tmp2(1),tmp2(end)];
xlim(x);

%for j = 1:2
subplot(3,1,1:2)
yyaxis left
y = yticks;
y = linspace(y(1),y(end),20);
%x = xlim;
for i = 1:size(allData.trialMeta.drugs,1)
    x(1) = allData.trialMeta.drugs{i,2};
    plot(x,[y(end-i+1),y(end-i+1)],'-k','linewidth',4)
    text(x(1),y(end-i+1),allData.trialMeta.drugs{i,1},'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',20)
end
%end
set(h1,'FontSize',20)
set(h2,'FontSize',20)
linkaxes([h1,h2],'x')
xlim([tmp2(1),tmp2(end)])
shg

%% define blocks to trial average
i = 0;
T = {};

while input('again?: ')
    i      = i+1;
    T{i,1} = input('name: ','s');
    T{i,2} = input('time start: ');
    T{i,3} = input('time end: ');
    T{i,4} = {};
    T{i,5} = {};
end

%% trial 1average
output_channel = 3;
rows = 1;
fr = 1e3;

t_vec = cellfun(@(s)(s.datetime),allData.trialData);
r_vec = cellfun(@(s)(s.inR),allData.trialData);

t               = allData.trialData{1}.time;
start_log       = diff(allData.trialData{1}.output(:,output_channel)) > 0;
stim_start      = seconds(t(start_log));
end_log         = diff(allData.trialData{1}.output(:,output_channel)) < 0;
stim_end        = seconds(t(end_log));


figure

for i = 1:size(T,1)
    idx =  t_vec > T{i,2} & t_vec < T{i,3} ;
    idx2 = find(diff(allData.trialData{1}.output(:,output_channel)) > 0,1);
    T{i,4} = cell2mat(cellfun(@(s)(s.scaledOutput),allData.trialData(idx),'UniformOutput',false)');
    T{i,4} = medfilt1(T{i,4},fr,'truncate'); %median filter the data
    T{i,5} = mean(T{i,4}(1:idx2,:),1); %(T{i,4}(end_log,:) - T{i,4}(start_log,:)) ./ (mean(T{i,4}(1:idx2,:),1)); %find rough conductance (Vf - Vi) / (Vm - Ve) = dg
    T{i,4} = T{i,4} - mean(T{i,4}(1:idx2,:),1); %baseline subtract the data
    
    tmp = range(T{i,4}) > 20;
    T{i,4}(:,tmp) = [];
    T{i,5}(:,tmp) = [];
    %T{i,4} = T{i,4} - mean(T{i,4}(1:idx2,:),1);
    
    v       = mean(T{i,4},2);
    sem     = std(T{i,4},[],2)/sqrt(size(T{i,4},2));
    
    h(i) = subplot(rows,ceil(size(T,1)/rows),i); hold on
    title(T{i,1})
    plot(t,v,'Color',[0,0,0],'linewidth',2)
    patch(seconds([t;flipud(t)]), [(v + sem); flipud(v - sem)],[1,0,0],'edgecolor','none','facealpha',0.5)
    plot(t,v,'Color',[0,0,0])
end

linkaxes(h)

y = ylim;
x = xlim;
if abs(y(2)) > abs(y(1))
    y_loc = y(2) - 0.05*(y(2) - y(1));
    va = 'top';
else
    y_loc = y(1) + 0.05*(y(2) - y(1));
    va = 'bottom';
end

for i = 1:size(T,1)

    patch(h(i),[stim_start,stim_start,stim_end,stim_end],[y(1),y(2),y(2),y(1)], [.1, 0, 1], 'edgecolor','none', 'facealpha',.1)
    text(h(i), x(2) - 0.05*(x(2)-x(1)), y_loc,sprintf('N = %i\n V_m = %i',size(T{i,4},2), round(mean(T{i,5}))),'HorizontalAlignment','right','VerticalAlignment',va,'FontSize',20) 
    ylim(y)
    xlim(x)
end
subplot(rows,ceil(size(T,1)/rows),1)
ylabel('voltage (mV))')


%% Plot Iont
plot_traces(allData,3); %inputs are data, and the channel (column) of output which corresponds to stimulus delivery
plot_effect(allData,3);

