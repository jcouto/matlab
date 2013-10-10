function [outvar] = process_dIV_file(ii,tt,plotvar)
%[outvar] = process_dIV_file(ii,tt,plotvar)
% ii - file number (from sorted list)  - defaul:t last file.
% tt - [tmin,tmax] interval used to compute the dIV - default: all noisy
% trace.
% plotvar - if defined produces a plot - default:yes


MAX_dIv_VALUE = 10;
MAX_dI_VOLT = -20;
MIN_dI_VOLT = -100;
% WINDOW = 100;
SAVE_FILE = 'dIV';
SAVE_FIG = 'dIV';

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

mask = spike_mask(Vall,diff(tall(1:2)));
Vrest = median(Vall(~mask))%(tall>0 & tall<=tprot(1))))
%Use only a stretch
t = tall(tall>tt(1)&tall<=tt(2));
V = Vall(tall>tt(1)&tall<=tt(2));
I = Iall(tall>tt(1)&tall<=tt(2));

dt = t(2)-t(1);
dVdt = (diff(V)*1e-3)./dt;

if ~isempty(plotvar)
    fig = figure(1);clf;
    ax(1) = axes('position',[.1,.1,.15,.25]);
    ax(2) = axes('position',[.35,.1,.25,.25]);
    ax(3) = axes('position',[.1,.45,.55,.25]);
    ax(4) = axes('position',[.1,.72,.55,.25]);
    ax(7) = axes('position',[.83,.15,.15,.2]); % isi distrib
    ax(5) = axes('position',[.85,.5,.15,.15]); % spike shape
    ax(6) = axes('position',[.75,.72,.2,.2]); % dVdt(V)
    ax(8) = axes('position',[.75,.43,.2,.22]); % Vm distrib
    ax(9) = axes('position',[.6,.15,.2,.2]); % fit result
    %     ax(4) = axes();
    axes(ax(1))
    
end

[C,capacitance.Ce, capacitance.var_I, p_C] = estimate_capacitance_from_noisy_trace(t,V,I,Vrest,plotvar);


if ~isempty(plotvar)
    axes(ax(2));
    xlabel('V (mV)');
    ylabel('I_m (pA)');
end

[eif, reif] = fit_rEIF_to_traces(t, V, I, C);

Im = (I - C*[dVdt(1),dVdt]);

% [dIV.raw_points.x, dIV.raw_points.y, Im, dI_V, dI_mu, dI_s, post_spk_data] = ...
%     extract_dIV(t, V, I, C, WINDOW, plotvar);
% 
% idx = (dI_V>MIN_dI_VOLT & dI_V<=MAX_dI_VOLT & -dI_mu/C <= MAX_dIv_VALUE );%
% % [x, f, resnorm, r,fit_output] = fit_EIF_to_dIV(dI_V(idx), dI_mu(idx),C);
% [x, f, resnorm, r,fit_output] = fit_rEIF_to_dIV(dI_V(idx), dI_mu(idx),...
%     post_spk_data.V,post_spk_data.Im,post_spk_data.t,C);

capacitance.C = C;
capacitance.Vrest = Vrest;

% Prepare output structure
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    expName = {'unknown'};
end
outvar.expName = expName{end};
outvar.capacitance = capacitance;
outvar.eif = eif;
outvar.reif = reif;
% outvar.dIV.Id = dI_mu;
% outvar.dIV.Ids = dI_s;
% outvar.dIV.v = dI_V;
% outvar.dIV.fit_idx = idx;
% outvar.eLIFparam = x;
% outvar.eLIFfun = f;
% /save(sprintf('dIV_%s.mat',expName{end}), '-struct', 'outvar');

if ~isempty(plotvar)
    axes(ax(9))
    plot(eif.dIV_V, -eif.dIV_Im/C,'k','linewidth',1)
    hold on
    plot(eif.dIV_V,eif.f(eif.param,eif.dIV_V),'r','linewidth',1)
   
    xlabel('V (mV)');
    ylabel('F(V) (mV/ms)');

    axes(ax(3));
    plot(t,V,'k')
    axis tight
    hold on
    Vm = integrate_eLIF(eif.param,t,I,C,V(1));
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
    
    bar(edge,n./trapz(edge,n),'facecolor','k','edgecolor','k')
    xlabel('V (mV)')
    %% Spike shapes
    axes(ax(5)),hold on
    [spk,spkw,tspkw, ~] = extract_spikes( V, [], t, 2, 5, 3);
    dspkw = diff(spkw, 1, 2)./(info.dt*1e3);
    plot(tspkw,spkw,'k','linewidth',0.5)
    mspkw = mean(spkw);
    plot(tspkw,mspkw,'r','linewidth',1.2)
    plot([tspkw(1),tspkw(end)],[0,0],'--','color',[.5,.5,.5])
    plot([tspkw(1),tspkw(end)],[-65,-65],'--','color',[.5,.5,.5])
    axis tight
    
    axes(ax(6)),hold on
    plot(spkw(:,1:end-1)',dspkw','k','linewidth',0.5)
    mdspkw = diff(mspkw)./(info.dt*1e3);
    plot(mspkw(1:end-1),mdspkw,'r','linewidth',1.2)
    axis tight
    xlabel('V (mV)')
    ylabel('dVdt (mV/ms)')
    axes(ax(7)),hold on
    edge = [1:5:300];
    [n] = histc(diff(spk*1e3),edge);
    bar(edge,n./trapz(edge,n),'facecolor','k','edgecolor','k')
    xlabel('isi (ms)')
    mean_isi = 1000/(length(spk)/(t(end)-t(1)));
    axis tight
    plot(mean_isi*[1,1],ylim(),'r')
    set(ax,'box','off','tickdir','out','linewidth',1,'color','none','ticklength',[0.03,0.03])
    set(ax(1),'xscale','log','yscale','log')
    set(ax(2),'xlim',[-80,-30],'ylim',[-3000,3000])
    set(ax(4),'xtick',[],'xcolor','w','ylim',[-3e3,3e3])
    set(ax(3:4),'xlim',[t(1),t(1)+1],'ticklength',[0.015,0.015])
    set(ax(5),'visible','off')
    linkaxes(ax(3:4),'x')
    set(fig,'paperposition',[0,0,18,15],'papersize',[18,15],'paperunits','centimeters')
    print(fig,'-dpdf',sprintf('dIV_%s.pdf',expName{end}))
    caption = sprintf(['Experiment %s. Measured noise statistics (Mean and ',...
        'standard deviation): pA'],expName{end})
    
    
end

function [t, V, I, metadata, info] = i_load_compensated_voltage(file,kfiles)
% Internal function to load voltage trace
V = [];
I = [];
metadata = [];
[ent, info] = load_h5_trace(file.path);
idx = find(strcmp('RealNeuron',{ent.name}));
idx = [idx, find(strcmp('AnalogInput',{ent.name}))];
if isempty(idx)
    idx = find(strcmp('ConductanceBasedNeuron',{ent.name}));
end
V = [ent(idx(1)).data];
t = linspace(0,info.tend,length(V));
% Do we need AEC?
idx = find(strcmp('Waveform',{ent.name}));
metadata = ent(idx).metadata;
% if there was a holding potential include it on the AEC current
idx = [idx,find(strcmp('Constant',{ent.name}))];
I = sum(vertcat(ent(idx).data),1);
if ~isempty(kfiles)
    % Is there AEC?    
    [~,k]  = min((file.date) - [kfiles.date]);
    Ke=load(kfiles(k).path);
    V = AECoffline(V,I,Ke);
end
