function thr = find_spike_threshold(x,stdmin)
% FIND_SPIKE_THRESHOLD finds the threshold for an extracellular trace
%
%   thr = FIND_SPIKE_THRESHOLD(x,stdmin)
%    Uses the formula from R. Quiroga, Z. Nadasdy, and Y. Ben-Shaul:
%       thr = stdmin*sigma_n ,with
%       sigma_n = median(|x|/0.6745)
% NOTE: 
%   Default stdmin is 4.
%   In the WaveClus batch scripts stdmin is set to 5.
% Adapted from WaveClus - Quiroga 2008
if ~exist('stdmin','var');stdmin = 5;end
thr = stdmin * median(abs(x))/0.6745;