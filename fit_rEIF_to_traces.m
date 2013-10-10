function [eif, reif] = fit_rEIF_to_dIV(t, V, I, C, plotvar)
% Fits a dynamic IV curve to a refractory exponential integrate and fire model
% neuron.
% It might be necessary to select the values of Y.
% Takes the voltage (X), mean transmembrane current (Y) and the Capacitance (C) as
% parameters
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
    plotvar = [1];
end
WINDOW = 200; % in ms 
MAX_FV = 10;

dt = diff(t(1:2));


[X, Y] = extract_dIV(t, V, I, C, WINDOW, plotvar);
% [x, f, resnorm, r,o] = fit_EIF_to_dIV(X(~isnan(Y)), Y(~isnan(Y)), C);
[x, f] = fit_EIF_to_dIV(X(~isnan(Y) & -Y/C<MAX_FV),...
    Y(~isnan(Y) & -Y/C<MAX_FV), C);

eif.dIV_V = X;
eif.dIV_Im = Y;
eif.param = x;
eif.f = f;

[~,spk_wave,tspk_wave, spkidx] = extract_spikes( V, [], t, 1, 3, 3);


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

NWINDOWS = cumsum(ceil(spk_width*3) + ...
    [0,5,5,10,30,30,30,30,50,50,100,100,...
    100,200,200,500,700,1000,1000]);

spkmask = spike_mask(V,dt);

dVdt = (diff(V)*1.0e-3)./dt;
Im = (I - C*[dVdt(1),dVdt]);
npoints = length(V);
edges = linspace(min(V),max(V),100);
refractory_param = nan(length(NWINDOWS)-1,4);
refractory_dIV_V = cell(length(NWINDOWS)-1,1);
refractory_dIV_Im = cell(length(NWINDOWS)-1,1);
for ii = 1:length(NWINDOWS)-1
        mask = false(size(V));
        window = int64((NWINDOWS(ii)/1000/dt:NWINDOWS(ii+1)/1000/dt));
        for jj = 1:length(spkidx)
            if (max(spkidx(jj)+window) < npoints)
                mask(spkidx(jj)+window) = 1;
            end
        end
%         disp([length(find(~mask | spkmask)),length(V)])
        maskedV = V; maskedV(~mask | spkmask) = nan;
        [X, Y, ~] = bin_samples(V(~isnan(maskedV)), Im(~isnan(maskedV)), edges);
        refractory_dIV_V{ii} = X;
        refractory_dIV_Im{ii} = Y;
        [tmpx] = fit_EIF_to_dIV(X(~isnan(Y) & -Y/C<MAX_FV),...
            Y(~isnan(Y) & -Y/C<MAX_FV), C, x);
        refractory_param(ii,:) = tmpx;
end

reif.dIV_V = refractory_dIV_V;
reif.dIV_Im = refractory_dIV_Im;
reif.param = refractory_param;

%% Plotting...
if ~isempty(plotvar)
    % Plot the spike shapes and the squares.
    SPIKES_TO_PLOT = 7;
    TPRE = 2;
    TPOST = NWINDOWS(5);
    nfig = figure();
    ax(1) = axes('position',[0.1,0.5,0.8,0.3]);
    idx = randperm(size(spk_wave,1));
    [spk_wave] = extractTriggeredTraces(V,spkidx(idx(1:SPIKES_TO_PLOT)),...
        int32(TPRE./1000/dt),...
        int32(TPOST./1000/dt));
    tspk_wave = linspace(-TPRE,TPOST,size(spk_wave,2));
    plot(tspk_wave,spk_wave,'k','linewidth',0.7)
    
    hold on
    tmp = min(min(spk_wave));
    plot([-2,0],tmp*[1,1],'k','linewidth',1.1)
    text(-1,tmp-0.5,'2ms','verticalalignment','bottom','horizontalalignment','center')
    axis tight;set(ax(1),'visible','off')
    % plot boxes
    for ii = 1:4
        rectangle('position',[NWINDOWS(ii),min(get(gca,'ylim')),diff(NWINDOWS([ii,ii+1]))-1,diff(get(gca,'ylim'))*.4],'linestyle',':')
    end
    
    for ii = 1:4
        ax(ii+1) = axes('position',[0.2+(ii-1)*.2,.8,.15,.15]);hold on
        annotation('textbox',[0.2+(ii-1)*.2,.9,.15,.05],...
            'string',sprintf('%d - %d ms',NWINDOWS(ii),NWINDOWS(ii+1)),...
            'edgecolor','none');
        
        plot(eif.dIV_V,eif.f(eif.param,eif.dIV_V),'k--')
        plot(reif.dIV_V{ii},-reif.dIV_Im{ii}/C,'ko','markersize',4)
        axis tight; ylim([min(ylim),10])
        
        plot(eif.dIV_V,eif.f(reif.param(ii,:),eif.dIV_V),'r-')
        plot(xlim,[0,0],'k:')
        set(ax,'visible','off')
        linkaxes(ax(2:end))
    end
end
% %%
% for ii = 1:4
%     ax(end+1) = axes('position',[0.1+(ii-1)*.22,.1,.17,.25]);hold on
%     xlabel('post-spike time (ms)')
% end
% %%
% 
% axes(ax(end-3))
% errorbar(W(1:end-1),1./tau_m,tau_m__,'ko')
% baseline = 1./x(1);
% plot(xlim,[1,1]*baseline,'b:')
% % sp = [0.2,10];
% s = fitoptions('Method','NonlinearLeastSquares');
% ff = fittype('a0 * (exp( -x / t0 ))','options',s);
% [params,G] = fit(W(1:end-1)',1./tau_m' - baseline,ff);
% plot(W(1:end-1),feval(params,(W(1:end-1)))+baseline,'r--')
% 
% % params = fit_decays(W(1:end-1),1./tau_m,lb,ub,sp,baseline);
% % plot(W(1:end-1),feval(params,(W(1:end-1)))+baseline,'b-')
% %
% ylabel('1/\tau_m (1/ms)')
% axes(ax(end-2))
% errorbar(W(2:end),E_m,E_m__,'ko')
% baseline = x(2);
% plot(xlim,[1,1]*baseline,'b:')
% s = fitoptions('Method','NonlinearLeastSquares');
% ff = fittype('a0 * (exp( -x / t0 ))','options',s);
% params = fit(W(2:end)',E_m' - baseline,ff);
% plot(W(1:end-1),feval(params,(W(1:end-1)))+baseline,'r--')
% ylabel('E_m (mV)')
% 
% axes(ax(end-1))
% errorbar(W(2:end),V_T,V_T__,'ko')
% baseline = x(2)+x(3);
% plot(xlim,[1,1]*baseline,'b:')
% 
% params = fit(W(2:end)',V_T'-baseline,ff)
% plot(W(2:end),feval(params,(W(2:end)))+baseline,'r--')
% ylabel('V_T (mV)')
% axes(ax(end))
% errorbar(W(2:end),delta_T,delta_T__,'ko')
% baseline = x(3);
% plot(xlim,[1,1]*baseline,':')
% lb = [-10,1e-2,-10,1e-2];
% ub = [30,1,30,1];
% sp = [6,0.01,0,0.01];
% params = fit(W(2:end)',delta_T'-baseline,ff);
% plot(W(2:end),feval(params,(W(2:end)))+baseline,'r--')
% 
% ylabel('\Delta_T (mV)')
%
% %%
% function params = fit_decays(X,Y,lb,ub,sp,baseline)
% % fits the decays with a double exponential
% 
% s = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',lb,...
%     'Upper',ub,...
%     'Startpoint',sp);
% ff = fittype('a0 * exp( -x / t0 ) + a1 * exp( -x / t1 )','options',s);
% [params,G] = fit(X',Y' - baseline,ff);
% 


