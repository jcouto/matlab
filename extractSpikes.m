function [spk,spk_w,t_spk_w] = extractSpikes(data, threshold, t, tpre, tpost, tdead, method)
% EXTRACTSPIKES Extracts the spikes from voltage data.
%
% [SPK, SPK_W, T_SPK_W] = extractSpikes( DATA, THRESHOLD, SRATE, TPRE, TPOST, TDEAD)
% Channels should be in collums.
%
%  Units:
%
% * DATA               - mV                            []
% * THRESHOLD          - mV                            [-10 mV]
% * TIME VECTOR|SRATE  - s | s-1                       [15e3 s-1]
% * TPRE               - ms (time window before spk)   [5ms]
% * TPOST              - ms (time window before spk)   [5ms]
% * TDEAD              - ms (min time between spk)     [2ms]
%
%  Extracts:
% * SPK (Timestamps of the action potentials)
% * SPK_W (Waveforms of the action potentials)%
% * T_SPK_W (Timevector for the waveforms of the action potentials)
% Note: SPK and SPK_W are cells if more than one row is provided to
%   DATA.

if ~exist('threshold','var'),threshold = -10;end
if ~exist('t','var')
    dt = 1./15e3;
    t  = (0:size(data,2)).*dt;
else
    if length(t)>1
        dt = t(2)-t(1);
    else
        dt = 1./t;
        t  = (0:size(data,1)).*dt;
    end
end
if ~exist('tpre','var'),tpre = 5;end
if ~exist('tpost','var'),tpost = 5;end
if ~exist('tdead','var'),tdead = 2;end
% if ~exist('method','var'),method = 'peak';end
% kink_threshold = 30; %mV/ms
if iscolumn(data)
    data = data';
end
N=size(data,1);

spk     = cell(N,1);
spk_w   = cell(N,1);
t_spk_w = cell(N,1);

parfor ii = 1:N
    idx         = argfindpeaks(data(ii,:), threshold, int32(tdead./1000/dt));
    spk_w{ii}   = extractTriggeredTraces(data(ii,:),idx,int32(tpre./1000/dt),int32(tpost./1000/dt));
    spk{ii}     = t(idx);
%     if strcmp(method,'kink')
        % Spikes are considered by the moment of crossing a mV ms-1 threshold
        % (20mV/ms is default)
%         for jj = 1:length(idx)
%             d_spk = diff(spk_w{ii}(jj,:))/(dt*1.e3);
%             idx(jj) = idx(jj) - abs(find(d_spk >= kink_threshold,1,'first') - ...
%                 argfindpeaks(spk_w{ii}(jj,:), threshold, tdead./1000/dt));
           
%         end
        % Get spikes triggered by threshold
%         spk_w{ii}   = extractTriggeredTraces(data(ii,:),idx,tpre./1000/dt,tpost./1000/dt);
%     end
    
end
% build time vector of spike waveforms
idx = find(~cellfun(@(x)isempty(x),spk_w),1,'first');
if ~isempty(idx) % but only if there are any spikes.
    t_spk_w = (t(1:size(spk_w{idx},2))-t(1))*1e3-tpre;
end

if (N < 2) % then return an array instead of a cell array.
    spk = spk{1};
    spk_w = spk_w{1};
end