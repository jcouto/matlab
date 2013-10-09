function [x, f, resnorm, r,o] = fit_rEIF_to_dIV(X, Y, postX, postY, postt, C, plotvar)
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

[x, f, resnorm, r,o] = fit_EIF_to_dIV(X(~isnan(Y)), Y(~isnan(Y)), C);



%%
if ~isempty(plotvar)
    % Plot the spike shapes and the squares before filtering second spikes.
    fig = gcf();
    nfig = figure();
    ax(1) = axes('position',[0.1,0.5,0.8,0.3]);
    idx = randperm(size(postX,1));
    plot(postt(postt<51),postX(idx(1:5),postt<51),'k','linewidth',0.7)
    
    hold on
    plot([-2,0],(min(min(postX(idx(1:5),postt<51))))*[1,1],'k','linewidth',1.1)
    text(-1,(min(min(postX(idx(1:5),postt<51)))-.5),'2ms','verticalalignment','bottom','horizontalalignment','center')
    axis tight;set(ax(1),'visible','off')
    % plot boxes
    rectangle('position',[5,min(get(gca,'ylim')),5-1,diff(get(gca,'ylim'))*.4],'linestyle',':')
    rectangle('position',[10,min(get(gca,'ylim')),10-1,diff(get(gca,'ylim'))*.4],'linestyle',':')
    rectangle('position',[20,min(get(gca,'ylim')),10-1,diff(get(gca,'ylim'))*.4],'linestyle',':')
    rectangle('position',[30,min(get(gca,'ylim')),20-1,diff(get(gca,'ylim'))*.4],'linestyle',':')
end

for ii  = 1:size(postX,1)
    [~,~,~,spk] = extract_spikes( postX(ii,:),[],postt./1000);
    mask = spike_mask(postX(ii,:),diff(postt(1:2))./1000);
    
    postX(ii,mask) = nan;
%     
%     if length(spk)>1
%         %         plot(postX(ii,:)),hold on
%         %         plot(spk,postX(ii,spk),'ko')
%         tmp = ceil(2000*(postt(2)-postt(1)));
%         for k = spk(2:end)'
%             if (k+tmp)>size(postX,2)
%                 postX(ii,k-5:end) = NaN;
%             else
%                 postX(ii,k-5:k+tmp) = NaN;
%             end
%         end
%         %         plot(postX(ii,:),'r'),hold on
%         %     pause,clf
%     end
end

W = [5,10,20,30,50,100,150,250,350,max(postt)];
ref_param = [];


for ii = 1:length(W)-1
%      for jj = 1:size(postX,1)
        
        tmpX = postX(:,postt>W(ii) & postt<W(ii+1));
        tmpY = postY(:,postt>W(ii) & postt<W(ii+1));
        [ntmpX, ntmpY, ~] = bin_samples(tmpX, tmpY,X);%linspace(min(tmpX(:)),max(tmpX(:)),20));
        idx = ~isnan(ntmpY);
        if length(ntmpX(idx))>5
            [tmpx, tmpf, ~, ~,~] = fit_EIF_to_dIV(ntmpX(idx), ntmpY(idx), C, x);
            ref_param = [ref_param;[W(ii),tmpx]];
            
        end
%         jj
%      end
end
%% [tau_m, E_m, delta_T, V_T]
% mean values
tau_m = nan(1,length(W)-1);
E_m = nan(1,length(W)-1);
delta_T = nan(1,length(W)-1);
V_T = nan(1,length(W)-1);
% standard errors
tau_m__ = nan(1,length(W)-1);
E_m__ = nan(1,length(W)-1);
delta_T__ = nan(1,length(W)-1);
V_T__ = nan(1,length(W)-1);

for ii = 1:length(W)-1
    tau_m(ii) = mean(ref_param(ref_param(:,1)==W(ii),2));
    tau_m___(ii) = std(ref_param(ref_param(:,1)==W(ii),2))./sqrt(length(ref_param(ref_param(:,1)==W(ii),2)));
    E_m(ii) = mean(ref_param(ref_param(:,1)==W(ii),3));
    E_m__(ii) = std(ref_param(ref_param(:,1)==W(ii),3))./sqrt(length(ref_param(ref_param(:,1)==W(ii),3)));
    delta_T(ii) = mean(ref_param(ref_param(:,1)==W(ii),4));
    delta_T__(ii) = std(ref_param(ref_param(:,1)==W(ii),4))./sqrt(length(ref_param(ref_param(:,1)==W(ii),4)));
    V_T(ii) = mean(ref_param(ref_param(:,1)==W(ii),5));
    V_T__(ii) = std(ref_param(ref_param(:,1)==W(ii),5))./sqrt(length(ref_param(ref_param(:,1)==W(ii),5)));
end


%%
if ~isempty(plotvar)
    pW = [5,10,20,30,50]; % window for the plots..
    for ii = 1:length(pW)-1
        ax(ii+1) = axes('position',[0.2+(ii-1)*.2,.8,.15,.15]);hold on
        annotation('textbox',[0.2+(ii-1)*.2,.9,.15,.05],...
            'string',sprintf('%d - %d ms',pW(ii),pW(ii+1)),...
            'edgecolor','none');
        
        plot(X,f(x,X),'k--')
        
        tmpX = postX(:,postt>pW(ii) & postt<pW(ii+1));
        tmpY = postY(:,postt>pW(ii) & postt<pW(ii+1));
        [ntmpX, ntmpY, ~] = bin_samples(tmpX, tmpY, X);
        plot(ntmpX(1:2:end),-ntmpY(1:2:end)/C,'ko','markersize',4)
        idx = ~isnan(ntmpY);
        [tmpx, tmpf, ~, ~,~] = fit_EIF_to_dIV(ntmpX(idx), ntmpY(idx), C, x);
        plot(linspace(min(ntmpX(idx)),-10),tmpf(tmpx,linspace(min(ntmpX(idx)),-10)),'r-')
        axis tight; ylim([min(ylim),10])
        
        plot(xlim,[0,0],'k:')
        set(ax,'visible','off')
    end
    linkaxes(ax(2:end))
end
%%
for ii = 1:4
    ax(end+1) = axes('position',[0.1+(ii-1)*.22,.1,.17,.25]);hold on
    xlabel('post-spike time (ms)')
end
%%

axes(ax(end-3))
errorbar(W(1:end-1),1./tau_m,tau_m__,'ko')
baseline = 1./x(1);
plot(xlim,[1,1]*baseline,'b:')
% sp = [0.2,10];
s = fitoptions('Method','NonlinearLeastSquares');
ff = fittype('a0 * (exp( -x / t0 ))','options',s);
[params,G] = fit(W(1:end-1)',1./tau_m' - baseline,ff);
plot(W(1:end-1),feval(params,(W(1:end-1)))+baseline,'r--')

% params = fit_decays(W(1:end-1),1./tau_m,lb,ub,sp,baseline);
% plot(W(1:end-1),feval(params,(W(1:end-1)))+baseline,'b-')
%
ylabel('1/\tau_m (1/ms)')
axes(ax(end-2))
errorbar(W(2:end),E_m,E_m__,'ko')
baseline = x(2);
plot(xlim,[1,1]*baseline,'b:')
s = fitoptions('Method','NonlinearLeastSquares');
ff = fittype('a0 * (exp( -x / t0 ))','options',s);
params = fit(W(2:end)',E_m' - baseline,ff);
plot(W(1:end-1),feval(params,(W(1:end-1)))+baseline,'r--')
ylabel('E_m (mV)')

axes(ax(end-1))
errorbar(W(2:end),V_T,V_T__,'ko')
baseline = x(2)+x(3);
plot(xlim,[1,1]*baseline,'b:')

params = fit(W(2:end)',V_T'-baseline,ff)
plot(W(2:end),feval(params,(W(2:end)))+baseline,'r--')
ylabel('V_T (mV)')
axes(ax(end))
errorbar(W(2:end),delta_T,delta_T__,'ko')
baseline = x(3);
plot(xlim,[1,1]*baseline,':')
lb = [-10,1e-2,-10,1e-2];
ub = [30,1,30,1];
sp = [6,0.01,0,0.01];
params = fit(W(2:end)',delta_T'-baseline,ff);
plot(W(2:end),feval(params,(W(2:end)))+baseline,'r--')

ylabel('\Delta_T (mV)')

%%

%%
function params = fit_decays(X,Y,lb,ub,sp,baseline)
% fits the decays with a double exponential

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,...
    'Upper',ub,...
    'Startpoint',sp);
ff = fittype('a0 * exp( -x / t0 ) + a1 * exp( -x / t1 )','options',s);
[params,G] = fit(X',Y' - baseline,ff);



