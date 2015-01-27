function [tau] = process_tau_folder(folder)
%[vi,fi,vm] = process_steps_folder(folder)
% This script computes the membrane time constant from a folder containing
% recorded traces.
% All files must have the same length.
% Not reading HOLDING value!
% Joao Couto, October 2013

if ~exist('plotvar','var');plotvar = 1;end
if ~exist('folder','var');folder = pwd;end

folder = cd(cd(folder));


% Parameters
% If TAU_TIME is not defined, uses the maximum off the trace
%TAU_TIME = 50; % in milliseconds.

%% Load all traces in folder
fprintf(1,'Loading traces from folder %s.\n',folder)
[files,kfiles] = list_h5_files(folder);
N = length(files);
if N < 1
    vi = [];
    fi = [];
    vm = [];
    fprintf(1,'Folder does not contain experiment files.\n');
    return
end
metadata = cell(N,1);
[t,V,I,metadata{1},info] = i_load_compensated_voltage(files(1),kfiles);
V = repmat(V,N,1);
I = repmat(I,N,1);
for ii = 2:length(files)
    [t,V(ii,:),I(ii,:),metadata{ii},info] = i_load_compensated_voltage(files(ii),kfiles);
end
% Discard traces with spikes

tmp_idx = find(max(V,[],2) > -20);
V(tmp_idx,:) = [];
I(tmp_idx,:) = [];
metadata(tmp_idx,:) = [];
% Search in the protocol until where to remove offset
t_prot = [0; cumsum(metadata{1,1}(:,1))];
tmp_idx = [find(t<=t_prot(1),1,'last'),find(t<=t_prot(2),1,'last')-1];
% Voltage trace with no offset and returns the membrane potential.
[Vc, o] = removeTraceOffset(V,tmp_idx);

dt = t(2)-t(1);
Vm = mean(Vc);
% Select the points to use in the fitting
idx0 = find(t>=t_prot(3),1,'first');

if ~exist('TAU_TIME','var')
    [~, TAU_TIME] = max(Vm(idx0:end));
else
    TAU_TIME = TAU_TIME/1000/dt;
end
idx1 = TAU_TIME + idx0;

srate = 1./dt;
% Single exponential fit
singleExp = fitExp(Vm(idx0:idx1), srate, 1);
ciSingle = confint(singleExp.F);
ciSingle = ciSingle(:,2);
% Double exponential fit
doubleExp = fitExp(Vm(idx0:idx1), srate, 2);
ciDouble = confint(doubleExp.F);
ciDouble = ciDouble(:,3:4);
taus = [doubleExp.F.b0,doubleExp.F.b1];
[taus,ndx] = sort(taus);
ciDouble = ciDouble(:,ndx);
% Select taus
if taus(1) < 15e-3 || (sum(isnan(ciSingle))>1)
    fprintf(1,'Using two exponentials.\n');
    ci = ciDouble;
else
    fprintf(1,'Using single exponential.\n');
    taus = singleExp.F.b0;
    ci = ciSingle;
end
% Save parameters to file.
holding_current = nan;

expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    try
        % Try to get the name from two folders above.
        [foldername,expName] = fileparts(fileparts(fileparts(pwd)));
    catch
        expName = {'unknown'};
    end
end
if iscell(expName)
   expName = expName{end};
end
%tmp = strsplit(folder,expName);
tmp = regexp(folder,expName,'split');
tmp = tmp{end};
tmp(tmp =='/') = '_';

Vm_offset = o;
tau.tau = taus;
tau.confidance = ci;
tau.Vm_offset = o;
tau.expName = expName;
tau.holding_current = holding_current;


appendix = sprintf('%s%s',expName,tmp);
DATANAME = 'tau.mat';

dataName = sprintf('%s/%s_%s',folder,appendix,DATANAME);
save(dataName,'-struct','tau')

% Plotting (run setFigureDefault before this.)
fig = figure(1);clf

% Plot the fit
ax = axes();
t0 = t(idx0);
T = t(idx0:idx1) - t0;
V = Vm(idx0:idx1);
hold on
plot(t(idx0:idx1)*1e3, Vc(:,idx0:idx1), 'Color', [.6,.6,.6]);
plot(t(idx0:idx1)*1e3, Vm(idx0:idx1), 'Color', [.3,.3,.3], 'LineWidth', 4);
fig_text = sprintf('%s\n\\tau_0 = %.2f ms',expName, taus(1)*1e3);
if length(taus) == 1
    plot((T+t0)*1e3, singleExp.F(T)+min(V), 'r', 'LineWidth', 2);
else
    plot((T+t0)*1e3, doubleExp.F(T)+min(V), 'r', 'LineWidth', 2);
    fig_text = sprintf('%s\n\\tau_1 = %.2f ms',fig_text, taus(2)*1e3);
end
text(min(xlim)+diff(xlim)./2, min(V)+1, fig_text)
axis tight;
axis([xlim, ylim+[-3,3]]);
xlabel('Time (ms)');
ylabel('Voltage (mV)');

% Plot the mean trace.
ax(2) = axes();
hold on
vmlims = [max(Vm),min(Vm)];
patch([t(idx0),t(idx0),t(idx1),t(idx1)]*1000,...
    [vmlims(2),vmlims(1),vmlims(1),vmlims(2)],[.7,.7,.7],...
    'edgecolor','none')
plot(t*1000,Vm,'Color', 'k', 'LineWidth', 1)
axis tight;
axis([xlim, ylim+[-1,1]]);
plot(min(xlim)+[.1,.1],min(Vm)+[0,2],'Color','k', 'LineWidth', 1)
plot(min(xlim)+[.1,30.1],min(Vm)+[0,0],'Color','k', 'LineWidth', 1)
text(min(xlim),min(Vm),'30ms','verticalalignment','top',...
    'horizontalalignment','left')
text(min(xlim),min(Vm),'2mV','verticalalignment','bottom',...
    'horizontalalignment','right')

% Plot the offset drift.
ax(3) = axes();
plot([1,length(o)],mean(o).*[1,1],'r',...
    'linewidth',2,'color',[.8,.3,.3])
hold on
plot(o,'ko','markersize',3.5,'markerfacecolor',[.4,.4,.4])
axis tight;
axis([xlim, ylim+[-5,5]]);
xlabel('Trace number');
ylabel('Initial voltage (mV)');

set(ax(1),'position',[0.1,0.1,0.5,0.6])
set(ax(2),'position',[0.1,0.72,0.8,0.25],'box','off','visible','off')
set(ax(3),'position',[0.65,0.1,0.25,0.6],'box','off','yaxislocation','right')

% Print figure and caption
set(gcf, 'Color', [1,1,1], 'PaperUnits', 'Inch', 'PaperPosition', [0,0,7,4],'PaperSize', [7,4]);

figName = sprintf('tau_%s.pdf',expName);
caption = sprintf(['Experiment %s - Short pulse protocol to compute the ',...
    'membrane time constant. Intensity of the pulse was %3.0f pA and the holding',...
    ' current %4.1f pA. Top - Average $V_m$ computed from %d traces.',...
    ' Bottom right - Stability, $V_m$ fluctuation across trials.',...
    ' Bottom left - Result of the fit. The offset was subtracted.'],...
    expName, min(I(:))-max(I(:)), holding_current, length(o));
if length(taus)==1
    caption = sprintf(['%s Used a single exponential to fit the data ',...
        'the time constant is %3.3fms.'],caption,taus(1)*1.e3);
else
    caption = sprintf(['%s Used two exponentials to fit the data ',...
        'the time constants are %3.3fms and %3.3fms.'],caption,taus(1)*1.e3,taus(2)*1.e3);
end

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
Vidx = idx(1);
V = [ent(Vidx).data];
t = linspace(0,info.tend,length(V));
% Do we need AEC?
idx = find(strcmp('Waveform',{ent.name}));
metadata = ent(idx).metadata;
% if there was a holding potential include it on the AEC current
idx = [idx,find(strcmp('Constant',{ent.name}))];
I = sum(vertcat(ent(idx).data),1);
if ~isempty(kfiles) & isempty(ent(Vidx).metadata)
    % Is there AEC?
    [~,k]  = min((file.date) - [kfiles.date]);
    Ke=load(kfiles(k).path);
    V = AECoffline(V,I,Ke);
end

