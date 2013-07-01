function [r0,r1,phi,conf,r,relspks,model_equation ] = compute_sinusoidal_modulation ...
    (spks, F, nBins, totalDuration, nPeriods, fitMeanRate, doPlot)
% [r0,r1,phi,conf] = extractSinusoidalModulation
%                        (spks, F, nBins, totalDuration, nPeriods, fitMeanRate, doPlot)
%
% Parameters:
%          spks - the spike times in seconds.
%             F - the frequency of the modulating current.
%         nBins - the number of bins to use for the histogram.
% totalDuration - the total duration of the stimulation.
%      nPeriods - the number of periods over which the fitting of the histogram
%                 of the spike times should be performed. Default value is 1.
%   fitMeanRate - whether to fit also the mean spiking rate. Default value is yes.
%        doPlot - whether to print intermediate results. Default value is no.
% 
% Returns:
%            r0 - the mean firing rate
%            r1 - the amplitude of the modulation of the firing rate.
%           phi - the phase of the modulation.
%          conf - the confidence intervals of the previous estimates.
%             r - firing rate histogram.
%       relspks - spikes relative to the sinusoid.
%         model - the model used in the fit.


% 
% Adapted from Daniele Linaro's scripts - June 2013.

if ~ exist('nPeriods','var')
    nPeriods = 1;
end
if ~ exist('fitMeanRate','var')
    fitMeanRate = 1;
end
relspks = mod(spks, nPeriods/F);
edges = linspace(0, nPeriods/F, nBins+1);
n = histc(relspks, edges);
r = n / (diff(edges(1:2)) * F * totalDuration);
r = r(:);
x = cumsum(diff(edges)) - 0.5*diff(edges(1:2));
model_equation = @(r0,r1,f,t,phi)(r0+r1*sin(2*pi*f*t+phi));

if fitMeanRate
    g = fittype('r0+r1*sin(2*pi*f*t+phi)','coeff',{'r0','r1','phi'},'problem','f','indep','t');
    model = fit(x(:),r(1:end-1),g,'StartPoint',[mean(r),(max(r)-min(r))/2,0],'problem',F);
    r0 = model.r0;
else
    r0 = mean(r(1:end-1));
    g = fittype('r0+r1*sin(2*pi*f*t+phi)','coeff',{'r1','phi'},'problem',{'f','r0'},'indep','t');    
    model = fit(x(:),r(1:end-1),g,'StartPoint',[(max(r)-min(r))/2,0],'problem',{F,r0});
end

r1 = model.r1;
phi = model.phi;
ci = confint(model);
conf = struct([]);
if fitMeanRate
    lbls = {'r0','r1','phi'};
else
    lbls = {'r1','phi'};
end
for k=1:length(lbls)
    conf(1).(lbls{k}) = ci(:,k)';
end

if r1 < 0
    r1 = abs(r1) + 1;
    phi = phi-pi;
    % fit the data again with the ``corrected'' values of r1 and phi as
    % starting points to obtain also the confidence intervals.
    if fitMeanRate
        model = fit(x(:),r(1:end-1),g,'StartPoint',[r0,r1,phi],'problem',F);
        r0 = model.r0;
    else
        model = fit(x(:),r(1:end-1),g,'StartPoint',[(max(r)-min(r))/2,0],'problem',{F,r0});
    end
    r1 = model.r1;
    phi = model.phi;
    ci = confint(model);
    for k=1:length(lbls)
        conf(1).(lbls{k}) = ci(:,k)';
    end
end

if exist('doPlot','var')
    x = linspace(0,edges(end),100);
    hold on;
    hndl = bar(edges, r, 0.8, 'histc');
    set(hndl,'FaceColor',[0.6,0.6,1]);
    plot(xlim,r0+[0,0],'g--','LineWidth',4);
    plot(x, r0+r1*sin(2*pi*F*x), 'm--', 'LineWidth', 4);
    plot(x, r0+r1*sin(2*pi*F*x+phi), 'k', 'LineWidth', 2);
    axis([0,edges(end),ylim]);
    set(gca, 'XTick', linspace(x(1),x(end),3), 'XTickLabel', {'0','0.5','1'});
    xlabel('Phase');
    ylabel('Rate (spikes/s)');
    title(sprintf('F = %.0f (%d spikes)', F, length(spks)));
end
