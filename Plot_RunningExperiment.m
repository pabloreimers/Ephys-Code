%% set parameters
fr = 10e-3 * 2e4; %define window for median filtering. with daqrate = 2e4, this is how many frames for the desired filtering window in objective time.
fs = 20; %fontsize for figures
ew = 50e-3 * 2e4;% define window for finding the effect size. Right now this is 50ms

%% Load in data
fprintf('select experiment data\n')         %ask user for data and load
[f,p]       = uigetfile('.mat','Select Data', 'C:\Users\ReimersPabloAlejandr\Documents\Code\pablo\Data\');
load([p,f])         

%% store raw data (into struct of D, called D)
%tmp         = regexp(p,'\','split');        %get the experiment name
%exp_name    = tmp{end-2};

idx         = any(allData.trialData{1}.output,1);   %find which output channel is used in the experiment

D.time      = allData.trialData{1}.time;            %store the template time for each trial
D.output    = allData.trialData{1}.output(:,idx);   %store the template output for each trial
D.voltage   = cell2mat(cellfun(@(x)(x.voltage),allData.trialData,'UniformOutput',false)');
D.current   = cell2mat(cellfun(@(x)(x.current),allData.trialData,'UniformOutput',false)');
D.soutput   = cell2mat(cellfun(@(x)(x.scaledOutput),allData.trialData,'UniformOutput',false)');
D.datetime  = cellfun(@(x)(x.datetime),allData.trialData,'UniformOutput',false)';
tmp         = cellfun(@(x)(x.datetime - allData.trialData{1}.datetime),allData.trialData,'UniformOutput',false)';
D.dt        = vertcat(tmp{:});
D.inR       = cell2mat(cellfun(@(x)(x.inR),allData.trialData,'UniformOutput',false)');
D.n         = size(allData.trialData,1);
D.drugs     = allData.trialMeta.drugs;

start_idx   = find(diff(D.output) > 0);
end_idx     = find(diff(D.output) < 0);
start_time  = D.time(start_idx);
end_time    = D.time(end_idx);

%% Process data
v   = medfilt1(D.soutput,fr,[],1,'truncate');
sem = std(v,[],2) / sqrt(size(v,2));
b   = mean(v(1:start_idx,:),1);

%% Plot mean scaled output and raster
figure(1); clf
subplot(3,1,1); hold on; title(exp_name);
patch([start_time,start_time,end_time,end_time],[min(v-b,[],'all'),max(v-b,[],'all'),max(v-b,[],'all'),min(v-b,[],'all')],[.1,0,1],'edgecolor','none','facealpha',.1)
plot(D.time,v - b,'color',[0,0,0,1/D.n])
patch([D.time;flipud(D.time)],[mean(v-b,2) + sem ; flipud(mean(v-b,2) - sem)],[1,0,0],'edgecolor','none', 'facealpha',0.5)
plot(D.time,mean(v-b,2),'r')
ylabel('V_m')
axis tight
%fontsize(gca,fs,'points')

subplot(3,1,2:3)
imagesc((v-b)')
colormap(bluewhitered)
xticks([])
ylabel('Trial #')

x = xlim;
for j = 1:size(D.drugs)
    [~,idx] = min(cellfun(@(x)(abs(x - D.drugs{j,2})),D.datetime));
    text(start_idx,idx,[D.drugs{j,1},'\rightarrow'],'HorizontalAlignment','right')
end
pos = get(gca,'Position');
c = colorbar;
ylabel(c,'V_m')
set(gca,'Position',pos)

%fontsize(gca,fs,'points')
%% Plot the effect size over trials
figure(2); clf
subplot(2,2,1) %plot the 
scatter(D.dt,v(end_idx,:) - b, 'filled','k')
ylabel('V_{end}')

subplot(2,2,2)
scatter(D.dt, max(v(start_idx:(end_idx+ew),:),[],1) - b,'filled','k')
ylabel('V_{max}')

%fit the rise and decay time
tau_r = zeros(D.n,1);
tau_d = zeros(D.n,1);


for i = 1:D.n
    tic
    [~,idx] = max(D.soutput(start_idx:(end_idx+ew),i)); %find the index of the max value in the effect size window
    x = seconds(D.time(start_idx+(1:idx)-1) - D.time(start_idx));
    y = -(v(start_idx+(1:idx)-1,i) - b(i));
    f = fit(x,y,'exp1');
    tau_r(i) = f.b;

    x = seconds(D.time((start_idx+idx):end) - D.time(start_idx+idx));
    y = v((start_idx+idx):end,i) - b(i);
    f = fit(x,y,'exp1');
    tau_d(i) = f.b;
    toc
end

subplot(2,2,3)
scatter(D.dt,tau_r,'filled','k'); ylabel('\tau_{rise}')

subplot(2,2,4)
scatter(D.dt,tau_d,'filled','k'); ylabel('\tau_{decay}')