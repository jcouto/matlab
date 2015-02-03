function [f,mm,vv,edges] = fano_factor(spikes,window)
% Computes fano factor from binary spike trains (N X time array)
%   - spikes: Array of binary spikes (N x time (1ms bins))
%   - window: to compute mean and variance 
%
% Returns:
%   -f: fano factor (spikes)
%   -mm: mean of spike counts (spikes)
%   -vv: variance of spike counts (spikes^2)

% First compute the spike cumulative sum

if ~exist('window','var')
    window = 100;
end
edges = 0:window:size(spikes,2)-1;

csum = cumsum(single(spikes), 2); % Compute cumulative spike count
% differentce between edges of window is the spike counts
counts = csum(:,edges(2:2:end) - edges(1) + 1) - csum(:,edges(1:2:end-1) - edges(1) + 1); 

mm = mean(counts);
vv = var(counts);
f = vv/mm;