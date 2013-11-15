function [onsets] = spike_onsets(v,dt,spks)
% [onsets] = spike_onsets(v,dt,spks)
% finds the onsets of the spikes for a given voltage trace (v) and
% spike indexes (spks) - optional. dt is the timestep of the recording (in s).
% Example: 
% [onsets] = spike_onsets(v,dt)
%
onsets = [];
if ~exist('v','var') || ~exist('dt','var');
    disp('Invalid number of inputs, usage: [onsets] = spike_onsets(v,dt);')
    return
end

if ~exist('spks','var'),[~,~,~,spkidx] = extract_spikes( v, [], dt);end

order = 2;
dvdt = diff(v,order)./((dt*1e3)^order);
if iscell(spkidx) % because it works only with single voltage traces...
    disp('Invalid number of traces or trace dimensions.')
    return
end
onsets = nan(size(spkidx));
previous = 1;
for ii = 1 : length(spkidx)
    % find the last peak of the voltage derivative
    previous = previous + ceil((spkidx(ii) - previous)*4/5);
    [vals,tmp] = findpeaks(dvdt(previous:spkidx(ii)-1),'minpeakheight',10);
    [~,idx] = max(vals);
    try
    onsets(ii) = previous + tmp(idx) - 1;
    catch
        fprintf(1,'spike_onsets: might of missed a peak here %d..\n',spkidx(ii))
    end
    previous = spkidx(ii);
end
onsets = onsets(~isnan(onsets));
