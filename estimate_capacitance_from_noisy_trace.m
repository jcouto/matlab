function [Cm,C,Variance_estimate,p] = estimate_capacitance_from_noisy_trace(t,V,I,Vmean,plotvar)
% Estimates the capacitance from an intracellular recording were a noisy
% current was injected.
% Accepts:
%   - t, time (s)    
%   - V, compensated membrane voltage (mV)
%   - I, injected current (pA)
%   - plotvar, arbitrary, if define plots creates a plot
%

%%% Parameters:
CMIN = 10;      %pF
CMAX = 1000;    %pF
CSTEP = 10;      %pF
DELTAV = 1;     %mV
p = [];
%%%
idx = find(I == I(1));
if ~exist('Vmean','var')
    Vmean = [];
end
if isempty(Vmean)
    Vmean = median(V(idx));
end
dt = t(2)-t(1);
dVdt = diff(V)./(dt*1.e3);

Ce = CMIN:CSTEP:CMAX;
Variance_estimate = nan(length(Ce),1);

idx = find(V(1:end-1)>=(Vmean-DELTAV) & V(1:end-1)<=(Vmean+DELTAV));

parfor ii = 1:length(Ce)
    Variance_estimate(ii) = var((I(idx)./Ce(ii)) - dVdt(idx));
end
[min_variance,C]=min(Variance_estimate);
Cm = Ce(C);

if exist('plotvar','var') && ~isempty(plotvar)
    
    p(1) = plot(Ce,Variance_estimate,'linewidth',0.7,'color','k');
    hold on;
    p(2) = plot(Cm,min_variance,'ko','markersize',4,'markerfacecolor',[.9,.3,.3]);
    xlabel('C_e (pF)')
    ylabel('Var[I_{in}/C_e-dV/dt]')
end
