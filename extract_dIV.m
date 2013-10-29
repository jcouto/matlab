function [dI_V, dI_mu, dI_s, X, Y, p] = extract_dIV(t, V, I, C, window, NBINS, plotvar)
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
if ~exist('NBINS','var')
    NBINS = 20;
end

p = [];
% Calculate Im(t)
dt = t(2)-t(1);
dVdt = (diff(V)*1.0e-3)./dt;

Im = (I - C*[dVdt(1),dVdt]);

% calculate dI(V) from trace without spikes
% Use threshold on dVdt
%[~,~,~, spkidx] = extract_spikes( dVdt, 20, t, 1, 3, 3);
% Use threshold on spike peaks
X = V(~isnan(V));
Y = Im(~isnan(V));

if (window ~= 0) 
    [~,~,~, spkidx] = extract_spikes( V, [], t, 1, 3, 3);
    
    window = ceil(window*1.0e-3/dt); %convert window to points
    tmpV = [V,nan(1,window)];    
    
    
    
    for jj =  1:length(spkidx)
        idx = spkidx(jj);
        tmpV(idx:idx+window) = nan;
    end
    
    tmpV = tmpV(1:end-window);
    X = tmpV(~isnan(tmpV));
    Y = Im(~isnan(tmpV));
end
% binning of the dI(V) curve
tmp = diff([min(V),max(V)]);
edges = linspace(min(V)+0.1*tmp,max(V)-0.3*tmp,NBINS);
% edges = linspace(min(V),max(V),NBINS);

[dI_V, dI_mu, dI_s] = bin_samples(X, Y, edges);

% %% Extract post spike parameters.
% 
% WPRE = int32(4*1e-3/dt);
% WPOST = int32(700*1e-3/dt);
% % remove fast spikes...
% % spk2remove = find(diff(spkidx)<WPOST);
% % for ii = flipud(spk2remove)'
% %     spkidx(ii-1) = [];
% %     spkidx(ii-1) = [];
% % end
% spkidx(find((spkidx-double(WPRE)) < 1)) = [];
% spkidx(find((spkidx+double(WPOST)) > length(V)-1)) = [];
% % if ((spkidx(1)-WPRE) < 1 ); spkidx(1) = [];end
% % if ((spkidx(end)+WPOST) > length(V)-1 ); spkidx(end) = [];end
% 
% spk_trig = nan(length(spkidx),WPRE+WPOST+1);
% Im_trig = nan(length(spkidx),WPRE+WPOST+1);
% tspk_trig = linspace(-double(WPRE)*1e3*dt,+double(WPOST)*dt*1e3,WPRE+WPOST+1);
% for jj =  1:length(spkidx)
%     idx = spkidx(jj);
%     spk_trig(jj,:)= V(idx-WPRE:idx+WPOST);
%     Im_trig(jj,:)= Im(idx-WPRE:idx+WPOST);
% end
% 
% post_spk.t = tspk_trig;
% post_spk.V = spk_trig;
% post_spk.Im = Im_trig;
% 
if exist('plotvar','var') && ~isempty(plotvar)
%     figure(fig)
    PLOT_POINTS = 0;
    ii = 0;
    if PLOT_POINTS
        ii = ii + 1;
        p(ii) = plot(X(1:10:end), Y(1:10:end),'o','markersize',3,...
            'markerfacecolor',[.5,.5,.5],'markeredgecolor',[.5,.5,.5]);
    end
    hold on
    ii = ii + 1;
    p(ii) = errorbar(dI_V,dI_mu,dI_s,'ok','markersize',3.5,...
        'markerfacecolor',[.9,.3,.3]);
    ii = ii + 1;
    p(ii) = plot(dI_V,dI_mu,'-r','linewidth',1);
    
end
% 
% % W = [5,10,20,30,50];
% % if exist('plotvar','var') && ~isempty(plotvar)
% %     fig = gcf;
% %     fig2 = figure();
% %     ax(1) = axes('position',[0.1,.6,.8,.3]);
% %     if size(spk_trig,1) > 5
% %         plot(tspk_trig,spk_trig(1:5,:),'k')
% %     else
% %         plot(tspk_trig,spk_trig,'k')
% %     end
% %     axis tight
% %     set(ax,'visible','off')
% %     % plot boxes
% %     rectangle('position',[5,min(get(gca,'ylim')),5-1,diff(get(gca,'ylim'))*.3],'linestyle','--')
% %     rectangle('position',[10,min(get(gca,'ylim')),10-1,diff(get(gca,'ylim'))*.3],'linestyle','--')
% %     rectangle('position',[20,min(get(gca,'ylim')),10-1,diff(get(gca,'ylim'))*.3],'linestyle','--')
% %     rectangle('position',[30,min(get(gca,'ylim')),20-1,diff(get(gca,'ylim'))*.3],'linestyle','--')
% %     
% %     ax(2) = axes('position',[0.43,.8,.1,.15]);hold on
% %     X = spk_trig(:,tspk_trig>5 & tspk_trig<10);
% %     Y = Im_trig(:,tspk_trig>5 & tspk_trig<10);
% %     [XX,YY] = bin_samples(X, Y, edges);
% %     [x, f, resnorm, r,fit_output] = fit_eLIF_to_dIV(X(~isnan(Y)), Y(~isnan(Y)),C);
% %     plot(XX,-YY/C,'r.');
% %     plot(XX,f(x,XX),'k-')
% %     ax(3) = axes('position',[0.55,.8,.1,.15]);hold on
% %     X = spk_trig(:,tspk_trig>10 & tspk_trig<20);
% %     Y = Im_trig(:,tspk_trig>10 & tspk_trig<20);
% %     [XX,YY] = bin_samples(X, Y, edges);
% %     [x, f, resnorm, r,fit_output] = fit_eLIF_to_dIV(X(:), Y(:),C);
% %     plot(XX,-YY/C,'r.');
% %     plot(XX,f(x,XX),'k-')    
% %     ax(4) = axes('position',[0.67,.8,.1,.15]);hold on
% %     X = spk_trig(:,tspk_trig>20 & tspk_trig<30);
% %     Y = Im_trig(:,tspk_trig>20 & tspk_trig<30);
% %     [XX,YY] = bin_samples(X, Y, edges);
% %     [x, f, resnorm, r,fit_output] = fit_eLIF_to_dIV(X(:), Y(:),C);
% %     plot(XX,-YY/C,'r.');
% %     plot(XX,f(x,XX),'k-')
% %     
% %     ax(5) = axes('position',[0.79,.8,.1,.15]);hold on
% %     X = spk_trig(:,tspk_trig>30 & tspk_trig<50);
% %     Y = Im_trig(:,tspk_trig>30 & tspk_trig<50);
% %     [XX,YY] = bin_samples(X, Y, edges);
% %     [x, f, resnorm, r,fit_output] = fit_eLIF_to_dIV(X(:), Y(:),C);
% %     plot(XX,-YY/C,'r.');
% %     %plot(XX,f(x,XX),'k-')
% %     set(ax(2:5),'yticklabel',{},'ticklength',[0.05,0.05],'box','on') 
% % end
% %%