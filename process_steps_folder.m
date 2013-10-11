function [vi,fi,vm] = process_steps_folder(folder)
%[vi,fi,vm] = process_steps_folder(folder)
% This function creates an fi and vi curves from a folder.
% "folder" input (default is the current working directory) is the folder name
% All files must have the same length.
% Joao Couto, October 2013

if ~exist('plotvar','var');plotvar = 1;end
if ~exist('folder','var');folder = pwd;end

folder = cd(cd(folder));

%% Load all traces in folder
fprintf(1,'Loading traces from folder %s.\n',folder)
[files,kfiles] = list_h5_files(folder);
N = length(files);
metadata = cell(N,1);
[t,V,I,metadata{1},info] = i_load_compensated_voltage(files(1),kfiles);
dt = info.dt;
V = repmat(V,N,1);
I = repmat(I,N,1);
for ii = 2:length(files)
    [t,V(ii,:),I(ii,:),metadata{ii},info] = i_load_compensated_voltage(files(ii),kfiles);
end

%% Extract spikes
[spks,spkw,tspkw] = extract_spikes(V,[],t);
isis = cellfun(@(x)diff(x),spks,'uniformoutput',0);

%% Find start and duration of stim
% NOTE: Preamble not supported.

tprot = cumsum(metadata{1}(:,1));
if length(tprot) > 3
    error('Protocol had more than 3 lines...')
end
tstim = tprot(1)+dt;
tpost = tprot(2);
dur = diff([tstim,tpost]);

%% Analize all steps:
%Select only traces with spikes during stim and no spikes before

vi.trials = false(N,1);
vi.v = nan(N,1);
vi.i = nan(N,1);
vi.v_sd = nan(N,1);
vi.v_sag = nan(N,1);
vi.dur = dur;
vi.rebound = nan(N,1);

fi.trials = false(N,1);
fi.f = nan(N,1);
fi.i = nan(N,1);
fi.n = nan(N,1);
fi.adapt = nan(N,1);
fi.dur = dur;
fi.fpost = nan(N,1);
fi.npost = nan(N,1);

vm.trials = false(N,1);
vm.v = nan(N,1);
vm.v_sd = nan(N,1);
vm.dur = tstim;

% Fraction of the stimulation to be used in the calculations of steady
% state
FRACTION = 3;

for ii = 1:N
    spk = spks{ii};
    isi = isis{ii};
    if isempty(spk(spk < tstim))
        vm.trials(ii) = 1;
        vm.v(ii) = mean(V(ii,t < tstim));
        vm.v_sd(ii) = std(V(ii,t < tstim));
        if ~isempty(spk((spk > tstim) & (spk < tpost)))
            fi.f(ii) = mean(1./isi((spk(1:end-1) > tpost - dur/FRACTION) & ...
                (spk(1:end-1) < tpost)));
            fi.i(ii) = mean(I(ii,(t > tpost-dur/FRACTION) & t < tpost)) -  ...
                mean(I(ii,t<tstim));
            fi.n(ii) = length(spk((spk > tstim) & (spk < tpost)));
            fi.fpost(ii) = mean(1./isi((spk(1:end-1) > tpost)));
            fi.npost(ii) = length(spk((spk(1:end) > tpost)));
            stim_isi = isi((spk(1:end-1) > tstim) & ...
                (spk(1:end-1) < tpost));
            fi.adapt(ii) = stim_isi(end)/stim_isi(1);
            fi.trials(ii) = 1;
        else
            vi.v(ii) = mean(V(ii,(t > tpost - dur/FRACTION) & ...
                (t < tpost)));
            vi.v_sd(ii) = std(((t > tpost - dur/FRACTION) & ...
                (t < tpost)));
            vi.v_sag(ii) = min(V(ii,(t > tstim ) & ...
                (t < tstim + dur/FRACTION)));
            vi.i(ii) = mean(I(ii,(t > tpost-dur/FRACTION) & t < tpost)) -  ...
                mean(I(ii,t<tstim));
            vi.trials(ii) = 1;
            
        end
    end
end

%% Fit fi and vi curves
fi.expr = 'a*x+b';
vi.expr = 'a*x+b';

[fi.coeff] = polyfit(fi.i(fi.trials),fi.f(fi.trials),1);
[vi.coeff] = polyfit(vi.i(vi.trials),(vi.v(vi.trials)-vm.v(vi.trials)),1);

%% Plots for fi, vi and vm
if exist('plotvar','var')
    cc = setFigureDefaults;
    fig = figure('visible','on');
    % Build frame
    ax(1) = axes('position',[0.1,0.75,.2,.2]);
    ylabel('Adaptation index') % adapt
    xlabel('Current (pA)')
    ax(2) = axes('position',[0.1,0.4,.2,.23]);
    ylabel('Frequency (Hz)')% fi
    xlabel('Current (pA)')
    ax(3) = axes('position',[0.1,0.1,.2,.23]);
    ylabel('\Delta Voltage')% vi
    xlabel('Current (pA)')
    ax(4) = axes('position',[0.4,0.5,.55,.45]);
    ylabel('Voltage (mV)')% fi and vi examples
    xlabel('Time (s)')
    ax(5) = axes('position',[0.4,0.1,.2,.25]);
    ylabel('Voltage (mV)')% Vm summary
    xlabel('Time (min)')
    ax(6) = axes('position',[0.7,0.1,.25,.25]);
    xlabel('Voltage (mV)')% Vm summary
    ylabel('dVdt (mV/ms)')
    
    label = 'ABCDEF';
    for ii = 1:length(ax)
        p = get(ax(ii),'position');
        ann(ii) = annotation('textbox',p+[-0.06,0.06,0,0],'string',label(ii));
    end
    set(ann,'verticalalignment','top',...
        'horizontalalignment','left',...
        'color','k','edgecolor','none',...
        'fontsize',10,'fontweight','bold')
    
    % Plot results
    if sum(fi.trials)
        axes(ax(1))
        plot(fi.i,fi.adapt,'ko','markerfacecolor',cc(1,:),...
            'markersize',2.5)
        ylim([0,nanmax(fi.adapt)*1.5])
        axes(ax(2))
        edges = [nanmin(fi.i)*0.9,nanmax(fi.i)*1.1];
        plot(edges,polyval(fi.coeff,edges),'k','linewidth',1)
        ylim([0,nanmax(fi.f)])
        plot(fi.i,fi.f,'ko','markerfacecolor',cc(1,:),...
            'markersize',2.5)
        plot(fi.i,fi.n./fi.dur,'k+')
        axis tight
    end
    if sum(vi.trials)
        axes(ax(3))
        edges = [nanmin(vi.i),nanmax(vi.i)];
        plot(edges,polyval(vi.coeff,edges),'k','linewidth',1)
        errorbar(vi.i,vi.v-vm.v,vi.v_sd,'ko',...
            'markerfacecolor',cc(1,:),'markersize',2.5)
        axis tight
        % Plot examples
        axes(ax(4)) % Pick 2 traces at random from fi and vi
        idx = find(vi.trials);
        idx = idx(randsample(length(idx),1));
        plot(t,V(idx,:),'k')
        axis tight
    end
    if sum(vm.trials)
        axes(ax(5)) % Plot the evolution of vm
        errorbar(([files.date]-min([files.date]))*60*24,vm.v,vm.v_sd,'ko',...
            'markerfacecolor',cc(1,:),'markersize',2.5)
        axis tight
        plot(xlim,nanmean(vm.v)+[0,0],'k--')
        ylim([-5,5]+nanmean(vm.v))
    end
    if sum(fi.trials)
        axes(ax(4))
        idx = find(fi.trials);
        idx = idx(randsample(length(idx),1));
        plot(t,V(idx,:),'color',cc(1,:))
        axis tight
        axes(ax(6)) % Plot the dv/V
        dspkw = diff(spkw{idx}, 1, 2)./(info.dt*1e3);
        idx2 = randsample(size(spkw{idx},1),10);
        plot(spkw{idx}(idx2,1:end-1)',dspkw(idx2,:)','color','k')
        plot(nanmean(spkw{idx}(:,1:end-1)),nanmean(dspkw(:,:)),'color','r')
        axis tight
        plot(xlim,[0,0],'k--')
        plot([0,0],ylim,'k--')
    end
end

%% Prepare caption
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    expName = {'unknown'};
end
expName = expName{end};

caption = sprintf(['Experiment %s. Performed %s. Computation of the cell FI',...
    ' and VI curves.'],expName,datestr(files(1).date));
if sum(fi.trials)
    caption = sprintf(['%s A - Adaptation index (last isi divided by the ',...
        'first) versus the injected current (pA) - mean $%1.2f\\pm%1.2f$.'],caption,...
        nanmean(fi.adapt),nanstd(fi.adapt));
    caption = sprintf(['%s B - Frequency-Current curve (N = %d traces). The red dots are the mean ISI ',...
        'in the last $%1.2f$s of the stimulus. The crosses the number of spikes divided ',...
        'by the duration. Firing frequency ranged from $%3.0f$ to $%3.0f$ Hz ($%3.0f$ to $%3.0f$ pA).',...
        ' The FI relation is $%3.1f$ Hz/nA.'],caption,sum(fi.trials),dur/FRACTION,nanmin(fi.f),nanmax(fi.f),...
        nanmin(fi.i),nanmax(fi.i),fi.coeff(1)*1e3);
else
    caption = sprintf(['%s There where no evoked spikes: A,B and F empty.'],caption);
end
if sum(vi.trials)
    caption = sprintf(['%s C - Voltage-Current curve (N = %d traces). The red dots are the mean voltage ',...
        'in the last $%1.2fs$ of the stimulus. The deflection ranged from $%3.0f$ to $%3.0f$ Hz',...
        ' ($%3.0f$ to $%3.0f$ pA). The membrane resistance is $%3.2f M\\Omega$.'],caption,...
        sum(vi.trials),dur/FRACTION,nanmin(vi.v),nanmax(vi.v),nanmin(vi.i),nanmax(vi.i),vi.coeff(1)*1e3);
else
    caption = sprintf(['%s There were no subthreshold trials so C is empty.'],caption);
end
if sum(vi.trials) | sum(fi.trials)
    caption = sprintf(['%s D - Example voltage traces of 2 examples.'],caption);
end
if sum(vm.trials) 
    caption = sprintf(['%s E - Stability of the membrane potential in the ',...
        'beginnign of each trial.'],caption);
end
if sum(fi.trials)
    caption = sprintf(['%s F - State plot of the first derivative versus the voltage during a spike.'],caption);
end
%% Save figures and data
tmp = strsplit(folder,expName);
tmp = tmp{end};
tmp(tmp =='/') = '_';
appendix = sprintf('%s%s',expName,tmp);

%%

FIGNAME = 'fi_vi.pdf';
DATANAME = 'fi_vi.mat';

dataName = sprintf('%s/%s_%s',folder,appendix,DATANAME);
save(dataName,'fi','vi','vm','files','kfiles','expName','folder')
figName = sprintf('%s/%s_%s',folder,appendix,FIGNAME);
%%
print(fig,'-dpdf',figName)
printFigWithCaption(figName,caption)
movefile([figName(1:end-4),'Caption.pdf'],figName)



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


