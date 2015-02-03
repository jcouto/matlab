function [bspks,interval] = binary_spiketrains(sp,interval,bin)
% [bspks] = binary_spiketrains(sp,interval,bin)
% Computes the spike counts from a spiketrain cell array (in seconds)
%   sp: array of spiketrains (works with 1d arrays only for the moment)
%   interval: the interval to binarize the spike trains (default is max of sp and min of sp)
%   bin: bin size (default 1ms)
%
%   bspks: are the binarized spiketrains (size of sp X time)
%

if ~exist('interval','var')
    interval = [];
end
if ~exist('bin','var')
    bin = 0.001;
end
%%%%%%%%%%Binarize spiketrains%%%%%%%%%%%%%%%
maxsp = nan(length(sp),1);
minsp = nan(length(sp),1);
for t = 1:length(sp)
    % because cellfun does not work when the trains are empty
    if ~isempty(sp{t})
        maxsp(t) = max(sp{t});
        minsp(t) = min(sp{t});
    end
end
maxsp = nanmax(maxsp);
minsp = nanmin(minsp);
if isempty(interval)
    interval = [minsp,maxsp];
end

edges = interval(1):bin:interval(2);
N  = length(edges);
M = length(sp);

bspks = zeros(M,N);
for m = 1:M
    bspks(m,:) = histc(sp{m},edges);
end
