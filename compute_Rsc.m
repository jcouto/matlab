function [rsc,r1,r2,zr1,zr2,bootstats] = compute_Rsc(sp1,sp2,window...
    ,method,nbootstrp)
% COMPUTE_RSC Computes Rsc according to Kohn et al. 2005
% [RSC,R1,R2,ZR1,ZR2,BOOTSTATS] = COMPUTE_RSC(SP1,SP2,WINDOW,METHOD,NBOOTSTR)
%
%   Computes Spike Count Correlations in a particular window (WINDOW) given
%the cellarrays that contain the spike times in each trial (SP1 and SP2). 
%   Discards the responses that are more than 3 times the standard deviation
% of either neurons mean response.
% METHOD can be 'rate' or 'zscore' (default). Specifies if the correlation
% is computed on the zscores or the spike counts themselves.
% NBOOTSTR is the number of bootstrap trials employed in the
% bootstrap (zero for no bootstrap).
% Outputs:
%   - RSC, the value of spike count correlation
%   - R1 and R1 the spikecounts in the defined window for the selected
%   trials.
%   - ZR1 and ZR2 the zscores of R1 and R2
%   - BOOTSTATS, the bootstrap results (distribution of the mean) only if
%   nbootstrp > 0
%

% either zscore or corr
if ~exist('method','var'); method = []; end
if isempty(method); method = 'zscore'; end
% Spike count window use all if not defined
if ~exist('window','var'); window = []; end
if isempty(window)
    window(1) = min(cellfun(@min,[sp1(:);sp2(:)]))-1;
    window(2) = max(cellfun(@max,[sp1(:);sp2(:)]))+1;
end

% Number of randomizations to compute significance using bootstrap method
% Zero for skipping randomization
if ~exist('nbootstrp','var'); nbootstrp = 0; end

% Get the spikes from the trials in the proper window.
r1 = cellfun(@(a)sum(a >= window(1) & a < window(2)),sp1);
r2 = cellfun(@(a)sum(a >= window(1) & a < window(2)),sp2);
filtvar = @(y)(y>=(mean(y)-3*std(y)) & y<=(mean(y)+3*std(y)));
idx = filtvar(r1) & filtvar(r2);
r1 =r1(idx);
r2 =r2(idx);

rsc = nan;
zr1 = [];
zr2 = [];
bootstats = [];
    
if length(r1) > 1
    zr1 = zscore(r1);
    zr2 = zscore(r2);
    switch method
        case 'zscore'
            x1 = zr1;
            x2 = zr2;
        case 'rate'
            x1 = r1;
            x2 = r2;
    end
    % Heavy lifting
    rsc = corr(x1,x2);
    if nbootstrp>0
        bootstats = bootstrp(nbootstrp,@(p,k)corr(p,k),x1,x2);
    end
end

