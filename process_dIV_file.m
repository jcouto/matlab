function [varargout] = process_dIV_file(ii,tt,plotvar)
% ii - file number (from sorted list)  - defaul:t last file.
% tt - [tmin,tmax] interval used to compute the dIV - default: all noisy
% trace.
% plotvar - if defined produces a plot - default:yes

[files,kfiles] = list_h5_files;

if ~exist('plotvar','var');plotvar = 1;end
if ~exist('ii','var');ii = length(files);end


[tall,Vall,Iall,metadata,info] = i_load_compensated_voltage(files(ii),kfiles);
% remove preamble or spontaneous.
tprot = unique(cumsum(metadata(:,1)));
if ~exist('tt','var')
    tt = [];
end
if isempty(tt)
    tt = [tprot(end) - tprot(end-1), tprot(end-1)];
end
Vrest = mean(Vall(tall>0 & tall<=tprot(1)));
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

[C,~,~,p_C] = estimate_capacitance_from_noisy_trace(t,V,I,Vrest,plotvar);

if ~isempty(plotvar),axes(ax(2));end

[rawVpoints, rawYpoints, Im, dI_V, dI_mu, dI_s, p_dIV] = ...
    extract_dIV(t, V, I, C, 50, plotvar);
idx = (dI_V>-100 & dI_V<=-40 & ~isnan(dI_mu)& -dI_mu/C <= 25 );%
[x, ~, f] = fit_eLIF_to_dIV(dI_V(idx), dI_mu(idx),C);


if ~isempty(plotvar)
    figure(2),clf
    plot(dI_V(idx), -dI_mu(idx)/C,'k')
    
    hold on
    plot(dI_V(idx),f(x,dI_V(idx)),'r')
    axes(ax(3));
    plot(t,V,'k')
    axis tight
    hold on
    Vm = integrate_eLIF(x,t,I,C,V(1));
    plot(t,Vm,'r')
    %     % eLIF response to DC current
    %     Vm2 = integrate_eLIF(x,t,150*ones(size(I)),C,V(1));
    %     plot(t,Vm2,'g')
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
    set(ax(3:4),'xlim',[t(1),t(1)+1])
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
    % if there was a holding potential include it on the AEC current
    idx = [idx,find(strcmp('Constant',{ent.name}))];
    I = sum(vertcat(ent(idx).data),1);
    [~,k]  = min((file.date) - [kfiles.date]);
    Ke=load(kfiles(k).path);
    V = AECoffline(V,I,Ke);
end
