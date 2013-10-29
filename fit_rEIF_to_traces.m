function [eif, reif,caption] = fit_rEIF_to_trace(t, V, I, C, plotvar)
% Fits a dynamic IV curve to a refractory exponential integrate and fire model
% neuron.
% It might be necessary to select the values of Y.
% Takes the voltage (X), mean transmembrane current (Y) and the Capacitance (C) as
% parameters
%   - plotvar use: [axes to plot dIV, figure to plot refractoriness]
%
% Returns:
%   - the result of the fit (parameters)
%   - resnorm the squared residuals
%   - f the function used to fit
% The result is organized in the following way: [tau_m, E_m, delta_T, V_T].
% For example:
%idx = (dI_V>-100 & dI_V<=-40 & ~isnan(dI_mu)& -dI_mu/C <= 25 );%
%[x] = fit_dIV_to_EIF(dI_V(idx), dI_mu(idx),C)
%%

if ~exist('plotvar','var')
    plotvar = [];
end
caption = '';
WINDOW = 200; % in ms 
MAX_FV = 15;
NBINS = 100;
USE_FIXED_BINS = 1; % Uses fixed bins after the spike (easier to average)...


dt = diff(t(1:2));

if length(plotvar)==2;axes(plotvar(1));end
[X, Y] = extract_dIV(t, V, I, C, WINDOW,NBINS, plotvar);
% [x, f, resnorm, r,o] = fit_EIF_to_dIV(X(~isnan(Y)), Y(~isnan(Y)), C);
[x, f] = fit_EIF_to_dIV(X(~isnan(Y) & -Y/C<MAX_FV),...
    Y(~isnan(Y) & -Y/C<MAX_FV), C);

eif.dIV_V = X;
eif.dIV_Im = Y;
eif.param = x;
eif.f = f;

[timestamps,spk_wave,tspk_wave, spkidx] = extract_spikes( V, [], t, 1, 3, 3);


%% Calculate the spike half width 
%(and set 3 times the spike half width as the starting point for the refractory)
spk_middle_volt = median(min(spk_wave,[],2) + (mean(max(spk_wave,[],2) - min(spk_wave,[],2))/2));
spk_width = 0;
for ii = 1:size(spk_wave,1) 
    tmp = tspk_wave(spk_wave(ii,:) > spk_middle_volt);
    spk_width = spk_width + diff([tmp(1),tmp(end)]);
end
spk_width = spk_width./size(spk_wave,1); % in miliseconds!
%%
% Compute the points to be calculated. And grab for each spike that is
% possible, the respective Im trace.
if USE_FIXED_BINS
    NWINDOWS = cumsum(5+[0,0,5,5,15,15,25,45,95,95]);
else % points depend on the spk_width
    NWINDOWS = cumsum(ceil(spk_width*3) + ...
        [2,5,5,10,30,30,30,30,50,50,100,100,...
        100,200,200,500,700,1000,1000]);
end
% Can only extract refractory properties when there are no spikes...
NWINDOWS = NWINDOWS(NWINDOWS < max(diff(timestamps*1e3))*.7);

spkmask = spike_mask(V,dt);

dVdt = (diff(V)*1.0e-3)./dt;
Im = (I - C*[dVdt(1),dVdt]);
npoints = length(V);
tmp = diff([min(V),max(V)]);
edges = linspace(min(V)+0.1*tmp,max(V)-0.4*tmp,NBINS);
refractory_param = nan(length(NWINDOWS)-1,4);
refractory_dIV_V = cell(length(NWINDOWS)-1,1);
refractory_dIV_Im = cell(length(NWINDOWS)-1,1);
for ii = 1:length(NWINDOWS)-1
        mask = false(size(V));
        window = int64((NWINDOWS(ii)/1000/dt:NWINDOWS(ii+1)/1000/dt));
        for jj = 1:length(spkidx)-1
            if (max(spkidx(jj)+window) < npoints) && (max(spkidx(jj)+window) < (1e-3 + spkidx((jj+1)))) 
                mask(spkidx(jj)+window) = 1;
            end
        end
%         disp([length(find(~mask | spkmask)),length(V)])
        maskedV = V; maskedV(~mask | spkmask) = nan;
        [X, Y, ~] = bin_samples(maskedV(~isnan(maskedV)), Im(~isnan(maskedV)), edges);
        refractory_dIV_V{ii} = X(3:end); % discard the first samples...
        refractory_dIV_Im{ii} = Y(3:end);
        [tmpx] = fit_EIF_to_dIV(X(~isnan(Y) & -Y/C<MAX_FV),...
            Y(~isnan(Y) & -Y/C<MAX_FV), C, x);
        refractory_param(ii,:) = tmpx;
end

reif.dIV_V = refractory_dIV_V;
reif.dIV_Im = refractory_dIV_Im;
reif.param = refractory_param;
reif.windows = NWINDOWS;

%% Plotting...
if ~isempty(plotvar)
    if length(plotvar)==2;figure(plotvar(2));end
    cc = setFigureDefaults;
    % Plot the spike shapes and the squares.
    SPIKES_TO_PLOT = 7;
    TPRE = 2;
    TPOST = NWINDOWS(5);
    ax(1) = axes('position',[0.1,0.1,0.8,0.5]);
    idx = randperm(size(spk_wave,1));
    [spk_wave] = extractTriggeredTraces(V,spkidx(idx(1:SPIKES_TO_PLOT)),...
        int32(TPRE./1000/dt),...
        int32(TPOST./1000/dt));
    tspk_wave = linspace(-TPRE,TPOST,size(spk_wave,2));
    plot(tspk_wave,spk_wave,'k','linewidth',0.7)
    
    hold on
    tmp = min(min(spk_wave));
    plot([-2,0],tmp*[1,1],'k','linewidth',1)
    text(-1,tmp-0.5,'2ms','verticalalignment','bottom','horizontalalignment','center')
    axis tight;set(ax(1),'visible','off')
    % plot boxes
    for ii = 1:4
        rectangle('position',[NWINDOWS(ii),min(get(gca,'ylim')),...
            diff(NWINDOWS([ii,ii+1]))-1,diff(get(gca,'ylim'))*.4],...
            'linestyle','--','edgecolor',cc(1,:))
    end
    
    for ii = 1:4
        ax(ii+1) = axes('position',[0.1+(ii-1)*.2,.6,.15,.30]);hold on
        annotation('textbox',[0.1+(ii-1)*.2,.85,.15,.05],...
            'string',sprintf('%d - %d ms',NWINDOWS(ii),NWINDOWS(ii+1)),...
            'edgecolor','none','horizontalalignment','left',...
            'verticalalignment','top');
        edges = linspace(min(V),max(V),NBINS*4);
        plot(edges,eif.f(eif.param,edges),'k--')
        plot(reif.dIV_V{ii},-reif.dIV_Im{ii}/C,'ko','markersize',3)
        tmp = find((-reif.dIV_Im{ii}/C)<15,1,'last');
        axis tight; ylim([min(ylim),15]),xlim([min(xlim),-30])%reif.dIV_V{ii}(tmp)])
        plot(edges,eif.f(reif.param(ii,:),edges),'r-')
        plot(xlim,[0,0],'k:')    
    end
    %linkaxes(ax(2:end))
    set(ax,'visible','off')
    caption = sprintf(['Example of %d spike waveforms and the obtained dynamic IV curves ',...
        'extracted in windows between %d to %d ms after each spike.'],...
        SPIKES_TO_PLOT,NWINDOWS(1),NWINDOWS(end));
end

