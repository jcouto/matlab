function [varargout] = process_dIV_file(ii,tt,plotvar)


[files,kfiles] = list_h5_files;

if ~exist('plotvar','var');plotvar = [];end
if ~exist('ii','var');ii = length(files);end


[tall,Vall,Iall,metadata,info] = i_load_compensated_voltage(files(ii),kfiles);
% remove preamble or spontaneous.
if ~exist('tt','var')
    tt = unique(cumsum(metadata(:,1)));
    tt = [tt(end) - tt(end-1), tt(end-1)];
elseif isempty(tt)
    tt = unique(cumsum(metadata(:,1)));
    tt = [tt(end) - tt(end-1), tt(end-1)];
end

%Use only a stretch
t = tall(tall>tt(1)&tall<=tt(2));
V = Vall(tall>tt(1)&tall<=tt(2));
I = Iall(tall>tt(1)&tall<=tt(2));
dt = t(2)-t(1);
dVdt = (diff(V)*1e-3)./dt;

if ~isempty(plotvar)
    fig = figure(1);clf;
    ax(1) = axes('position',[.1,.1,.2,.25]);
    ax(2) = axes('position',[.4,.1,.3,.25]);
    ax(3) = axes('position',[.1,.43,.6,.25]);
    ax(4) = axes('position',[.1,.72,.6,.25]);
    ax(7) = axes('position',[.75,.1,.2,.25]); % isi distrib
    ax(5) = axes('position',[.85,.2,.15,.15]); % spike shape
    ax(6) = axes('position',[.75,.72,.2,.22]); % dVdt(V)
    ax(8) = axes('position',[.75,.43,.2,.22]); % Vm distrib
    
%     ax(4) = axes();
    axes(ax(1))
    
end
[C,~,~,p_C] = estimate_capacitance_from_noisy_trace(t,V,I,-65,plotvar)

if ~isempty(plotvar),axes(ax(2));end
[X, Y, Im, dI_V, dI_mu, dI_s, p_dIV] = extract_dIV(t, V, I, C, 50, plotvar);
idx = (dI_V>-100 & dI_V<=-40 & ~isnan(dI_mu)& -dI_mu/C <= 25 );%
[x] = fit_eLIF_to_dIV(dI_V(idx), dI_mu(idx),C)

figure(2),clf
plot(dI_V(idx), -dI_mu(idx)/C,'k')
f = @(a,X)(1.0/a(1))*(a(2) - X + (a(3) * exp((X - a(4))/a(3))))
hold on
plot(dI_V(idx),f(x,dI_V(idx)),'r')
if ~isempty(plotvar)
    axes(ax(3));
    plot(t,V,'k')
    axis tight
    hold on
    axes(ax(4));
    plot(t,Im,'k')
    axis tight
    %% 
    axes(ax(8)),hold on
    edge = [min(V)-1:.5:max(V)+1];
    [n] = histc(V,edge);
    bar(edge,n,'facecolor','k','edgecolor','k')
    %% Spike shapes
    axes(ax(5)),hold on
    [spk,spkw,tspkw, ~] = extract_spikes( V, [], t, 2, 5, 3);
    dspkw = diff(spkw, 1, 2)./(info.dt*1e3);
    plot(tspkw,spkw,'k','linewidth',0.5)
    mspkw = mean(spkw);
    plot(tspkw,mspkw,'r','linewidth',1.2)
    plot([tspkw(1),tspkw(end)],[0,0],'--','color',[.5,.5,.5])
    xlim([tspkw(1),tspkw(end)])
    axes(ax(6)),hold on
    plot(spkw(:,1:end-1)',dspkw','k','linewidth',0.5)
    mdspkw = diff(mspkw)./(info.dt*1e3);
    plot(mspkw(1:end-1),mdspkw,'r','linewidth',1.2)
    axis tight
    axes(ax(7)),hold on
    [n] = histc(diff(spk*1e3),[1:5:300]);
    bar([1:5:300],n,'facecolor','k','edgecolor','k')
    
    mean_isi = 1000/(length(spk)/(t(end)-t(1)));
    axis tight
    plot(mean_isi*[1,1],ylim(),'r')
    set(ax,'box','off','tickdir','out','linewidth',1,'color','none')
    set(ax(1),'xscale','log','yscale','log')
    set(ax(2),'xlim',[-80,-30],'ylim',[-3000,3000])
    set(ax(4),'xtick',[],'xcolor','w','ylim',[-3e3,3e3])
    set(ax(3:4),'xlim',[t(1),t(1)+0.3])
    linkaxes(ax(3:4),'x')
end
    
function [t, V, I, metadata, info] = i_load_compensated_voltage(file,kfiles)
% Internal function to load voltage trace
V = [];
I = [];
metadata = [];
[ent, info] = load_h5_trace(file.path);
idx = find(strcmp('RealNeuron',{ent.name}));
idx = [idx, find(strcmp('AnalogInput',{ent.name}))];
V = [ent(idx(1)).data];
t = linspace(0,info.tend,length(V));
% Do we need AEC?
if isempty(ent(idx(1)).metadata) && ~isempty(kfiles)
    % Is there I?
    idx = find(strcmp('Waveform',{ent.name}));
    metadata = ent(idx).metadata;
    idx = [idx,find(strcmp('Constant',{ent.name}))];
    I = sum(vertcat(ent(idx).data),1);
    [~,k]  = min((file.date) - [kfiles.date]);
    Ke=load(kfiles(k).path);
    V = AECoffline(V,I,Ke);
end