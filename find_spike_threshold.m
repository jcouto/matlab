function thr = find_spike_threshold(x,stdmin)
% Finds the threshold for an extracellular trace
%    Uses the formula from R. Quiroga, Z. Nadasdy, and Y. Ben-Shaul:
%    thr = stdmin*sigma_n ,with
%    sigma_n = median(|x|/0.6745)
% NOTE: 
%   Default stdmin is 4.
%   In the WaveClus batch scripts stdmin is set to 5.
% Adapted from WaveClus - Quiroga 2008
if ~exist('stdmin','var');stdmin = 4;end
thr = stdmin * median(abs(x))/0.6745;