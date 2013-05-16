function [X, Y, Im, dI_V, dI_mu, dI_s, p] = extract_dIV(t, V, I, C, window,plotvar)
% Extracts the dynamic IV curve from raw data
% [X, Y, dI_V, dI_mu, dI_s, Im] = extract_dIV(t, V, I, C, window,plotvar)
% Takes as inputs:
%   t - time(s)
%   V - Vm the compensated membrane voltage
%   C - membrane capacitance estimate (use "estimate_capacitance_from_noisy_trace")
%   I - injected current (pA)
%   window - the window to remove after each spike (ms)
%   plotvar - if defined a plot is generated.
% The outputs are:
%   X - points of V without the time after the spikes given by "window"
%   Y - points of Im without the time after the spikes given by "window"
%   Im - the transmembrane current (pA)
%   dI_V - binned samples of X
%   dI_mu - binned samples of Y
%   dI_s - std of the distribution in dI_mu
%   p - the plot objects (requires plotvar defined)

if ~exist('C','var') 
    C = estimate_capacitance_from_noisy_trace(t,V,I);
elseif isempty(C)
    C = estimate_capacitance_from_noisy_trace(t,V,I);
end
if ~exist('window','var')
    window = 200;%ms
end
p = [];
% Calculate Im(t)
dt = t(2)-t(1);
dVdt = (diff(V)*1.0e-3)./dt;

Im = (I - C*[dVdt(1),dVdt]);

% calculate dI(V) from trace without spikes
% Use threshold on dVdt
[~,~,~, spkidx] = extract_spikes( dVdt, 20, t, 1, 3, 3);
% Use threshold on spike peaks
% [~,~,~, spkidx] = extract_spikes( V, [], t, 1, 3, 3);

window = ceil(window*1.0e-3/dt); %convert window to points

tmpV = [V,nan(1,window)];

for jj =  1:length(spkidx)
    idx = spkidx(jj);
    tmpV(idx:idx+window) = nan;
end

tmpV = tmpV(1:end-window);

% binning of the dI(V) curve
edges = linspace(-90,-30,100);

X = tmpV(~isnan(tmpV));
Y = Im(~isnan(tmpV));

[dI_V, dI_mu, dI_s] = bin_samples(X, Y, edges);


if exist('plotvar','var') && ~isempty(plotvar)
    p(1) = plot(X, Y,'o','markersize',3,...
        'markerfacecolor',[.5,.5,.5],'markeredgecolor',[.5,.5,.5]);
    hold on
    p(2) = errorbar(dI_V,dI_mu,dI_s,'ok','markersize',3.5,...
        'markerfacecolor',[.9,.3,.3]);
    p(3) = plot(dI_V,dI_mu,'-k','linewidth',1);
end