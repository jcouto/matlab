function [Cm,Ce,Variance_estimate,p] = estimate_capacitance_from_noisy_trace(t,V,I,Vrest,deltav,plotvar)
% Estimates the capacitance from an intracellular recording were a noisy
% current was injected.
% Accepts:
%   - t, time (s)    
%   - V, compensated membrane voltage (mV)
%   - I, injected current (pA)
%   - Vrest, the resting membrane potential.
%   - plotvar, arbitrary, if define plots creates a plot
%
% NOTE:
% The resolution is dictated by the variable CSTEP, which can be set in
% the fuction code (for ploting only).


%%% Parameters:
CMIN = 1;      %pF
CMAX = 1000;    %pF
CSTEP = 1;      %pF
DELTAV = 1;     %mV
CGUESS = 250.0;

p = [];
%%%
if ~exist('Vrest','var')
    Vrest = [];
end

if ~exist('deltav','var')
    deltav = DELTAV;
end

if isempty(Vrest)
    mask = spike_mask(V,diff(t(1:2)));
    Vrest = median(V(~mask));
end
dt = t(2)-t(1);
dVdt = diff(V)./(dt*1.e3);
idx = find(V(1:end-1)>=(Vrest-deltav) & V(1:end-1)<=(Vrest+deltav));
Fx = @(x)var((I(idx)./x) - dVdt(idx));
[Cm] = fminsearch(Fx,CGUESS);


if nargout>1 || exist('plotvar','var')
    Ce = CMIN:CSTEP:CMAX;
    Variance_estimate = nan(size(Ce));
    parfor ii = 1:length(Ce)
        Variance_estimate(ii) = Fx(Ce(ii));
    end
    [min_variance]=min(Variance_estimate);
end
if exist('plotvar','var') && ~isempty(plotvar)
    p(1) = plot(Ce,Variance_estimate,'linewidth',0.7,'color','k');
    hold on;
    p(2) = plot(Cm,min_variance,'ko','markersize',4,'markerfacecolor',[.9,.3,.3]);
    xlabel('C_e (pF)')
    ylabel('Var[I_{in}/C_e-dV/dt]')
end
