function [outvar] = process_eif_file(fname,tt,C,plotvar)
%[outvar] = process_eif_file(fname,tt,plotvar)
%     fname - name of the file to be analysed
%     tt - [tmin,tmax] interval used to compute the dIV - default: all noisy
%     trace.
%     plotvar - if defined produces a plot - default:yes


MAX_dIv_VALUE = 10;
MAX_dI_VOLT = -20;
MIN_dI_VOLT = -100;
DATANAME = 'eif.mat';


if ~exist('plotvar','var');plotvar = 1;end

%% Find the file and compensate if needed.

[foldername,fname] = fileparts(fname);

[files,kfiles] = list_h5_files;
fileidx = find(cellfun(@(x)~isempty(x),cellfun(@(x)(strfind(cell2mat(x),...
    cell2mat(regexp(fname,'(\d)','match')))),...
    regexp({files.path},'(\d)','match'),'uni',0)));
if isempty(fileidx);
    outvar = [];
    fprintf(1,'Could not find file: %s\n',fname);
    return
end
[foldername,filename] = fileparts(files(fileidx).path);
fprintf(1,'Processing file: %s\n',files(fileidx).basename);
foldername = cd(cd(foldername));
[tall,Vall,Iall,metadata,info] = i_load_compensated_voltage(files(fileidx),kfiles);

%% Trim data
% Remove preamble or spontaneous.
tprot = unique(cumsum(metadata(:,1)));
if ~exist('tt','var')
    tt = [];
end
if isempty(tt)
    tt = [tprot(end) - tprot(end-1), tprot(end-1)];
end

dt = tall(2)-tall(1);

mask = spike_mask(Vall, dt);
Vrest = median(Vall(~mask));%(tall>0 & tall<=tprot(1))))
% Use the first part of the trace to get Vrest
%if there are no spikes 

% if isempty(find(Vall(tall<tt(1))>-20))
%     Vrest = mean(Vall(tall<tt(1)));
% end


t = tall(tall>tt(1) & tall<=tt(2));
V = Vall(tall>tt(1) & tall<=tt(2));
I = Iall(tall>tt(1) & tall<=tt(2));
dVdt = (diff(V)*1e-3)./dt;

if ~isempty(plotvar)
    %%
    cc = setFigureDefaults;
    label = 'ABCDEF';
    
    fig(1) = figure(1);clf;clear tmp
    tmp(1) = axes('position',[.1,.7,.8,.2]); % Im trace
    xlabel('t (s)'); ylabel({'Transmembrane','current (I_m) (pA)'}); hold on
    tmp(2) = axes('position',[.1,.4,.8,.2]); % V trace
    
    xlabel('t (s)'); ylabel('Voltage (mV)'); hold on
    
    tmp(3) = axes('position',[.1,.1,.2,.2]); % Capacitance
    tmp(4) = axes('position',[.4,.1,.2,.2]); % Im vs V
    xlabel('V (mV)'); ylabel('F(V) (mV/ms)'); hold on
    tmp(5) = axes('position',[.7,.1,.2,.2]); % fit result
    xlabel('V (mV)');ylabel('I_m (pA)');
    ax{1} = tmp;
    
    fig(2) = figure(2);clf;clear tmp
    tmp(1) = axes('position',[.1,.55,.35,.35]);
    ylabel('P(V)'); xlabel('Voltage (mV)'); hold on
    tmp(2) = axes('position',[.55,.55,.35,.35]);
    ylabel('P(isi)'); xlabel('interspike interval (ms)'); hold on
    tmp(3) = axes('position',[.1,.1,.35,.35]);
    xlabel('time (ms)'); ylabel('Voltage (mV)'); hold on
    tmp(4) = axes('position',[.55,.1,.35,.35]); % spike shape
    xlabel('Voltage (mV)'); ylabel('dV/dt (mV/ms)'); hold on
    ax{2} = tmp;
    for jj = 1:length(fig)
        figure(fig(jj))
        for ii = 1:length(ax{jj})
            p = get(ax{jj}(ii),'position');
            ann{jj}(ii) = annotation('textbox',p+[-0.06,0.06,0,0],'string',label(ii));
        end
        set(ann{jj},'verticalalignment','top',...
            'horizontalalignment','left',...
            'color','k','edgecolor','none',...
            'fontsize',10,'fontweight','bold')
        
    end
    %%
    fig(3) = figure(3);clf;
    
    axes(ax{1}(3))
    caption = cell(3,1);
end

%% Estimate capacitance
deltav = 1;
[Cm,capacitance.Ce, capacitance.var_I] = estimate_capacitance_from_noisy_trace(t,V,I,Vrest,deltav,plotvar);

if ~exist('C','var');C = [];end
if isempty(C);C = Cm;end

%% Get eif parameters
[eif, reif,caption{3}] = fit_rEIF_to_traces(t, V, I, C,[ax{1}(4),fig(3)]);

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


%% Prepare output and save data
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    expName = {'unknown'};
end
expName = expName{end};

outvar.capacitance = capacitance;
outvar.eif = eif;
outvar.reif = reif;


tmp = regexp(foldername,expName,'split');
tmp = tmp{end};
tmp(tmp =='/') = '_';
appendix = sprintf('%s%s',filename,tmp);

save(sprintf('%s/%s_%s',foldername,appendix,DATANAME), '-struct', 'outvar');

%% Plot
if ~isempty(plotvar)
    axes(ax{1}(5))
    idx = find(-eif.dIV_Im/C < 40,1,'last');
    idx = 1:idx;
    plot(eif.dIV_V(idx), -eif.dIV_Im(idx)/C,'k','linewidth',1)
    edges = linspace(min(V),max(V),100);
    axis tight
    plot(edges,eif.f(eif.param,edges),'r','linewidth',1)
%     tmpidx1 = find(-eif.dIV_Im/C < 20 & eif.dIV_V < 0,1,'last');
    
%     ylim([min(-eif.dIV_Im/C),10])
%     xlim([eif.dIV_V([1]),-30])
%     axes(ax{1}(4))
%     xlim(eif.dIV_V(idx))
    ylim([min(ylim)-1,10])
    set(ax{1}(4),'xlim',get(ax{1}(5),'xlim'))
    
    axes(ax{1}(2));
    idx = (t>tprot(end-2)-0.05 & t<tprot(end-2)+3);    
    plot(t(idx),V(idx),'k')
    axis tight
    hold on
    Vm = integrate_EIF(eif.param,dt*1e3,I,C,V(1));
    plot(t(idx),Vm(idx),'r')
    %     % eLIF response to DC current
    %     Vm2 = integrate_eLIF(x,t,150*ones(size(I)),C,V(1));
    %     plot(t,Vm2,'g')
    axes(ax{1}(1));
    
    plot(t(idx),Im(idx),'k')
    axis tight
    
    set(ax{1},'box','off','tickdir','out','linewidth',0.7,'color','none')
    set(ax{1}(3),'xscale','log','yscale','log')
%     set(ax{1}(4),'ylim',[-3000,3000])
    set(ax{1}(1),'xtick',[],'xcolor','w','ylim',[min(Im)*0.07,max(Im)*0.17])
    
    linkaxes(ax{1}(1:2),'x')
    
    caption{1} = sprintf(['Experiment %s - %s: estimation of the exponential ',...
        'integrate and fire parameters. A - Transmembrane current (',...
        '$I_m(t) = I_{in} - CdV/dt$). B - Example of a recorded trace ',...
        '(compensated with AEC). C - Minimization of the $var[I_{in}/C_e-dV/dt]$',...
        ' to estimate the cell capacitance (%3.1fpF). Values of V at ',...
        '$%2.1f \\pm %2.1f mV$ where used. ',...
        'D - Im (pA) versus Vm (mV) to estimate the dynamic IV curve. ',...
        'E - Dynamic IV curve plotted as $F(V)=-I_d/C$ together',...
        ' with the exponential fit. EIF parameters:(%3.1f,%3.1f,%3.1f,%3.1f).',...
        ''],expName,filename,capacitance.C,Vrest,deltav,eif.param(1),...
        eif.param(2),eif.param(3),eif.param(4));
    
    %%
    [spk,spkw,tspkw, ~] = extract_spikes( V, [], t, 2, 5, 3);
    axes(ax{2}(1)),hold on
    edge = (min(V)-1:.5:max(V)+1);
    [n] = histc(V,edge);
    bar(edge,n./trapz(edge,n),'facecolor','k','edgecolor','k')
    xlabel('V (mV)')
    
    axes(ax{2}(2)),hold on
    edge = [1:5:300];
    isi = diff(spk*1e3);
    [n] = histc(isi,edge);
    
    bar(edge,n./trapz(edge,n),'facecolor','k','edgecolor','k')
    mean_isi = 1000/(length(spk)/(t(end)-t(1)));
    axis tight
    plot(mean_isi*[1,1],ylim(),'r')
    
    axes(ax{2}(3)),hold on
    dspkw = diff(spkw, 1, 2)./(info.dt*1e3);
    plot(tspkw,spkw,'k','linewidth',0.5)
    mspkw = mean(spkw);
    plot(tspkw,mspkw,'r','linewidth',1.2)
    plot([tspkw(1),tspkw(end)],[0,0],'--','color',[.5,.5,.5])
    plot([tspkw(1),tspkw(end)],[-65,-65],'--','color',[.5,.5,.5])
    axis tight
    axes(ax{2}(4)),hold on
    plot(spkw(:,1:end-1)',dspkw','k','linewidth',0.5)
    mdspkw = diff(mspkw)./(info.dt*1e3);
    plot(mspkw(1:end-1),mdspkw,'r','linewidth',1.2)
    axis tight
    xlabel('V (mV)')
    ylabel('dVdt (mV/ms)')
    
    
    caption{2} = sprintf(['Experiment %s - %s: Action potential and membrane voltage',...
        ' statistics. A - histogram of the membrane voltage ($%3.1f \\pm %3.1f mV$).',...
        ' B - Distribution of the interspike intervals (mean: $%3.3f$ ; median: ',...
        '$%3.1f ms$). C - Average spike shape. D - State plot of the first derivative',...
        ' of the voltage versus the voltage during a spike.'],expName,filename,nanmean(V),...
        nanstd(V),mean(isi),median(isi));
    caption{3} = sprintf(['Experiment %s - %s: Estimation of refractory properties. %s'],...
        expName,filename,caption{3});
    set(fig(1:2),'paperposition',[0,0,18,15],'papersize',[18,15],'paperunits','centimeters')
    set(fig(3),'paperposition',[0,0,18,10],'papersize',[18,10],'paperunits','centimeters')
    %% Save figures and data
    
    FIGNAME = {'eif.pdf','Vm_and_spike.pdf','refractory_prop.pdf'};
    
    for ii = 1:length(fig)
        figName = sprintf('%s/%s_%s',foldername,appendix,FIGNAME{ii});
        print(fig(ii),'-dpdf',figName)
        printFigWithCaption(figName,caption{ii})
        movefile([figName(1:end-4),'Caption.pdf'],figName)
    end
end

function [t, V, I, metadata, info] = i_load_compensated_voltage(file,kfiles)
% Internal function to load voltage trace

FILTERDATA = 0;
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
    if FILTERDATA
        V = filter_data(V,[],5000,1.0/info.dt);
    end
end
