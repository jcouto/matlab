function mask = spike_mask(v,dt,onsets,N)
% mask = spike_mask(v,dt,onsets)
% Returns a matrix of zeros that is one at the location of the spikes.

mask = [];
if ~exist('v','var') || ~exist('dt','var');
    disp('Invalid number of inputs, usage: [mask] = spike_mask(v,dt);')
    return
end
if ~exist('N','var')
    N = 5;% number of points (extra) to remove before and after the spikes
end

if ~exist('onsets','var'),[onsets] = spike_onsets( v, dt);end

mask = logical(zeros(size(v)));

dvdt = diff(v)./(dt);


for ii = 1:length(onsets)
    if ((ii+1)>=length(onsets))
        argnext = length(dvdt);
    else
        argnext = onsets(ii+1);
    end
    argnext = onsets(ii) + ceil(argnext-onsets(ii))/3;
    tmp = find(dvdt(int64(onsets(ii):argnext))<0,1,'first');
    tmp2 = find(dvdt(int64(onsets(ii)+tmp:argnext))>0,1,'first');
%     [vals,tmp] = findpeaks(-dvdt(int64(onsets(ii)):int64(argnext)));
%     [~,idx] = max(vals);
    mask(int64(onsets(ii)-N:onsets(ii)+tmp2+N)) = 1;
end
% figure('visible','on')
% t = linspace(0,60,length(v));
%  plot(t,v,'b')
%  hold on
% 
% vv = v;
% vv(~mask) = nan;
% plot(t,vv,'r')
% pause
% close(gcf)
%     

