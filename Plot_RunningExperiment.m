%% set parameters
fr = 30e-3 * 2e4; %define window for median filtering. with daqrate = 2e4, this is how many frames for the desired filtering window in objective time.
fs = 20; %fontsize for figures
ew = 200e-3 * 2e4;% define window for finding the effect size. Right now this is 50ms
sw = 0e-3 * 2e4;

tau_flag = false;
%% Load in data
fprintf('select experiment data\n')         %ask user for data and load
[f,p]       = uigetfile('.mat','Select Data', 'C:\Users\ReimersPabloAlejandr\Documents\Code\pablo\Data\');
load([p,f])         

%% store raw data (into struct of D, called D)
tmp         = regexp(p,'\','split');        %get the experiment name
exp_name    = tmp{end-2};

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
plot(D.time,v - b,'color',[0,0,0,10/D.n])
patch([D.time;flipud(D.time)],[mean(v-b,2) + sem ; flipud(mean(v-b,2) - sem)],[1,0,0],'edgecolor','none', 'facealpha',0.5)
plot(D.time,mean(v-b,2),'r')
ylabel('V_m')
axis tight
fontsize(gca,fs,'points')

subplot(3,1,2:3)
imagesc((v-b)')
colormap(bluewhitered)
xticks([])
ylabel('Trial #')

x = xlim;
for j = 1:size(D.drugs)
    [~,idx] = min(cellfun(@(x)(abs(x - D.drugs{j,2})),D.datetime));
    text(start_idx,idx-1,[D.drugs{j,1},'\rightarrow'],'HorizontalAlignment','right')
end
pos = get(gca,'Position');
c = colorbar;
ylabel(c,'V_m')
set(gca,'Position',pos)

fontsize(gca,fs,'points')

%% Fit rise and decay times
if tau_flag
tau_r = zeros(D.n,1);
tau_d = zeros(D.n,1);

ft = fittype('(a ./ (1+exp(-b*(x-c)))) + d','independent','x','dependent','y');
opts = fitoptions('Method','NonlinearLeastSquares','Display','Off');

psp  = v(start_idx:(end_idx+ew),:) - v(start_idx,:);     %extract the post synaptic potential
[m,idx] = max(abs(psp),[],1);       %find the max effect in the psp (baseline subtracted)
s       = sign(psp(sub2ind(size(psp),idx,1:D.n)));           %extract the sign of that effect

for i = 1:D.n
    tic
    opts.StartPoint = [s(i)*m(i),100,0.05, 0]; %initialize the sigmoidal fit where max value is carrying capacity(a), inflection point is end of stim(c), and imperically chosen 100 for scaling(b), and 0 offset (d).
    xq= start_idx+(1:idx(i))-1;
    x = seconds(D.time(xq) - D.time(start_idx));
    y = v(start_idx+(1:idx(i))-1,i) - b(i);
    f = fit(x,y,ft,opts);
    tau_r(i) = f.b;
    
    x = seconds(D.time((start_idx+idx(i)):end) - D.time(start_idx+idx(i)));
    y = v((start_idx+idx(i)):end,i) - b(i);
    f = fit(x,y,'exp1');
    tau_d(i) = f.b;
    toc
end
end

%% Plot the effect size over trials
psp     = v((start_idx+sw):(start_idx+ew),:);     %extract the post synaptic potential
psp1    = psp(1,:);

figure(2); clf
h(1) = subplot(2,3,1); 
scatter(D.dt,max(psp-psp1,[],1), 'filled','k')
ylabel('V_{max}')
drug_label(D)

h(2) = subplot(2,3,2);
scatter(D.dt, min(psp-psp1,[],1),'filled','k')
ylabel('V_{min}')
drug_label(D)

h(2) = subplot(2,3,3);
scatter(D.dt, mean(psp-psp1,1),'filled','k')
ylabel('V_{avg}')
drug_label(D)

h(5) = subplot(2,3,5);
scatter(D.dt,b,'filled','k'); ylabel('V_{rest}')
drug_label(D)

if tau_flag
h(4) = subplot(2,3,4);
tau_r(isoutlier(tau_r,'percentiles',[.5,99.5])) = nan;
scatter(D.dt,tau_r,'filled','k'); ylabel('\tau_{rise}')
drug_label(D)

h(5) = subplot(2,3,5);
tau_d(isoutlier(tau_d,'percentiles',[.5,99.5])) = nan;
scatter(D.dt,tau_d,'filled','k'); ylabel('\tau_{decay}')
drug_label(D)
end

h(6) = subplot(2,3,6);
D.inR(isoutlier(D.inR,'percentiles',[.5,99.5])) = nan;
scatter(D.dt,D.inR,'filled','k'); ylabel('\Omega_{input}')
drug_label(D)

linkaxes(h,'x')
%% Plot example fit for tau
if tau_flag
i = 100;
psp  = v(start_idx:(end_idx+ew),:) - v(start_idx,:);     %extract the post synaptic potential
[m,idx] = max(abs(psp),[],1);       %find the max effect in the psp (baseline subtracted)
s       = sign(psp(sub2ind(size(psp),idx,1:D.n)));           %extract the sign of that effect

opts.StartPoint = [s(i)*m(i),100,0.05, 0]; %initialize the sigmoidal fit where max value is carrying capacity(a), inflection point is end of stim(c), and imperically chosen 100 for scaling(b), and 0 offset (d).
xq= start_idx+(1:idx(i))-1;
x = seconds(D.time(xq) - D.time(start_idx));
y = v(start_idx+(1:idx(i))-1,i) - b(i);
f = fit(x,y,ft,opts);

figure(3); clf
%show rise time fit
subplot(1,2,1)
plot(D.time(xq),y)
hold on
plot(D.time(xq), (f.a ./ (1 + exp(-f.b * (x - f.c)))) + f.d)
tmp1 = xlim;
tmp2 = ylim;
text(tmp1(1),tmp2(2),sprintf(['$$f(x) = \\frac{a}{1 + e^{-b(x-c)}} - d$$\n' ...
    '                  a = %.2f\n' ...
    '                   b = %.2f\n' ...
    '                   c = %.2f\n' ...
    '                   d = %.2f\n'],...
    f.a,f.b,f.c,f.d),'HorizontalAlignment','left','verticalAlignment','top','Interpreter','Latex')

%show decay fit
x = seconds(D.time((start_idx+idx(i)):end) - D.time(start_idx+idx(i)));
y = v((start_idx+idx(i)):end,i) - b(i);
f = fit(x,y,'exp1');
ylabel(['V_m of trial ',num2str(i)])
fontsize(gca,fs,'points')

subplot(1,2,2)
plot(D.time((start_idx+idx(i)):end),y)
hold on
plot(D.time((start_idx+idx(i)):end),f.a.*exp(f.b.*x))
tmp1 = xlim;
tmp2 = ylim;
text(tmp1(1),tmp2(2),sprintf(['$$f(x) = ae^{bx}$$\n' ...
    '                  a = %.2f\n' ...
    '                   b = %.2f\n'],...
    f.a,f.b),'HorizontalAlignment','left','verticalAlignment','top','Interpreter','Latex')
legend('data','fit')
fontsize(gca,fs,'points')
end

%% show trial average by drug condition and individual traces
if ~strcmp('saline',D.drugs{1,1})
    D.drugs = [{'saline',D.datetime{1}};D.drugs];
end

names       = {size(D.drugs(:,1))};
start_times = {size(D.drugs(:,1))};
end_times   = {size(D.drugs(:,1))};

figure(4); clf; clear h h2
subplot(2,1,1)
hold on
for i = 1:length(D.inR)
    scatter(D.datetime{i},D.inR(i),'k','filled')
end
drug_label(D,0)
subplot(2,1,2)
hold on
tmp = mean(psp-b,1);
for i = 1:length(D.inR)
    scatter(D.datetime{i},tmp(i),'k','filled')
end
drug_label(D,0)

for i = 1:size(D.drugs,1)
    names{i} = D.drugs{i,1};

    a = input([D.drugs{i,1},' start time: ']);
    if isa(a,'datetime')
        start_times{i} = a;
    else
        start_times{i} = D.drugs{i,2};
    end

    a = input([D.drugs{i,1},' end time: ']);
    if isa(a,'datetime')
        end_times{i} = a;
    elseif i < size(D.drugs,1)
        end_times{i} = D.drugs{i+1,2};
    else
        end_times{i} = D.datetime{end};
    end
end

for i = 1:size(names,2)
    idx = cellfun(@(x)(x > start_times{i} && x < end_times{i}),D.datetime);
    
    m   = mean(v(:,idx)-b(idx),2);
    sem = std(v(:,idx),[],2)/sqrt(sum(idx));

    h(i) = subplot(2,size(D.drugs,1),i); hold on
    plot(D.time,m,'k','LineWidth',2)
    patch([D.time;flipud(D.time)],[m+sem;flipud(m-sem)],[1,0,0],'EdgeColor','none','FaceAlpha',.5)
    axis tight
    title(names{i})
    xticklabels([])

    h2(i) = subplot(2,size(D.drugs,1),i+size(D.drugs,1));
    [~,idx] = min(cellfun(@(x)(abs(end_times{i} - x)),D.datetime));
    plot(D.time,v(:,idx-1),'k','Linewidth',1)
    axis tight
end
linkaxes(h)
ylabel(h(1),'Trial Average Voltage (mV)')
y = ylim(h(1));
x = xlim(h(1));
ylabel(h2(1),'Singel Trial Voltage (mV)')

for i = 1:size(names,2)
    patch(h(i),[start_time,start_time,end_time,end_time],[y(1),y(2),y(2),y(1)],[.1,0,1],'edgecolor','none','facealpha',0.1)
    chi=get(h(i), 'Children');
    set(h(i),'Children',flipud(chi))
    idx = cellfun(@(x)(x > start_times{i} && x < end_times{i}),D.datetime);
    text(h(i), x(1), y(1),...
        sprintf('N = %i\n$$V_{rest}$$ = %.1fmV\n$$R_{input} = %.0fM\\Omega$$',sum(idx), mean(b(idx)),mean(D.inR(idx),'omitnan')),...
        'HorizontalAlignment','left','VerticalAlignment','top','FontSize',20,'Interpreter','latex') 
    fontsize(h(i),fs,'points')

    subplot(2,size(D.drugs,1),i+size(D.drugs,1))
    x_tmp = xlim;
    y_tmp = ylim;
    patch(h2(i),[start_time,start_time,end_time,end_time],[y_tmp(1),y_tmp(2),y_tmp(2),y_tmp(1)],[.1,0,1],'edgecolor','none','facealpha',0.1)
    chi=get(h2(i), 'Children');
    set(h2(i),'Children',flipud(chi))
    fontsize(h2(i),fs,'points')
end

%% Plot how principal components vary over time
[coeff,score,latent] = pca(v);
latent = latent/sum(latent); %normalize the variance so it sums to 1

n = sum(latent > .1); % find how many priniciple comopnents carry more than 5 percent of the variance in the dataset. they are probably not noise then
figure(5); clf
c = make_colors(n);
for i = 1:n
    subplot(2,n,i)
    plot(D.time,score(:,i),'LineWidth',2,'Color',c(i,:))
    title(sprintf('PCA %i, \\sigma = %i%%',i,round(100*latent(i))))

    subplot(2,n,i+n)
    scatter(D.dt,coeff(:,i),'filled','k')
    drug_label(D)
end


