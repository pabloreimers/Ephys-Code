function drug_label(D,start_time)
if ~exist('start_time','var')
    start_time = D.datetime{1};
end

for i = 1:size(D.drugs,1)
    xline(D.drugs{i,2} - start_time, ':')
    tmp = ylim;
    h = text(D.drugs{i,2} - start_time,tmp(2),D.drugs{i,1},'HorizontalAlignment','right','VerticalAlignment','top');
    set(h,'Rotation',90);
end

drawnow