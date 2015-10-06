function [spk,spk_w,t_spk_w,spk_idx] = extract_spikes(data, threshold, t, tpre, tpost, tdead, method)
% EXTRACT_SPIKES Extracts the spikes from voltage data.
%
% [spk,spk_w,t_spk_w,spk_idx] = extract_spikes(data, threshold, t, tpre, tpost, tdead, method)
% Channels should be in collums.
%
%  Units:
%
% * DATA               - mV                            []
% * THRESHOLD          - mV                            [-10 mV]
% * TIME VECTOR|DT     - s | s                       [1./30e3 s]
% * TPRE               - ms (time window before spk)   [5ms]
% * TPOST              - ms (time window before spk)   [5ms]
% * TDEAD              - ms (min time between spk)     [2ms]
%
% Note: if threshold is empty, the threshold is selected from the media of
% all points that are larger than 100 mV/ms.
%  Extracts:
% * SPK (Timestamps of the action potentials)
% * SPK_W (Waveforms of the action potentials)%
% * T_SPK_W (Timevector for the waveforms of the action potentials)
% * SPKI_DX (Index of the action potentials)
% Note: SPK and SPK_W are cells if more than one row is provided to
%   DATA.

if iscolumn(data)
    data = data';
end

if ~exist('threshold','var'),threshold = -10;end; tt = threshold;

if ~exist('t','var')
    dt = 1./30e3;
    t  = (0:size(data,2)).*dt;
else
    if length(t)>1
        dt = t(2)-t(1);
    else
        dt = t;
        t  = (0:size(data,2)).*dt;
    end
end
if ~exist('tpre','var'),tpre = 5;end
if ~exist('tpost','var'),tpost = 5;end
if ~exist('tdead','var'),tdead = 2;end
if ~exist('method','var'),method = 'peak';end
% kink_threshold = 30; %mV/ms

N=size(data,1);

spk     = cell(N,1);
spk_w   = cell(N,1);
t_spk_w = cell(N,1);
spk_idx     = cell(N,1);

for ii = 1:N
    if isempty(threshold)
        dVdt = diff(data(ii,:)*1e-3)./dt;
        tt = min(median(data(ii,dVdt>100)));
    end
    idx         = argfindpeaks(data(ii,:), tt, int32(tdead./1000/dt));
    
    tmp   = extractTriggeredTraces(data(ii,:),idx,int32(tpre./1000/dt),int32(tpost./1000/dt));
    spk_w{ii} = tmp;
    spk{ii}     = t(idx);
    spk_idx{ii} = idx;
end
% build time vector of spike waveforms
idx = find(~cellfun(@(x)isempty(x),spk_w),1,'first');
if ~isempty(idx) % but only if there are any spikes.
    t_spk_w = (t(1:size(spk_w{idx},2))-t(1))*1e3-tpre;
end

if (N < 2) % then return an array instead of a cell array.
    spk = spk{1};
    spk_w = spk_w{1};
    spk_idx = spk_idx{1};
end