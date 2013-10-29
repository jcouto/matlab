function [var] = process_eif_folder(foldername,C,FORCE,plotvar)

if ~exist('C','var');C = []; end
if ~exist('FORCE','var');FORCE = 0; end
if ~exist('plotvar','var');plotvar = 1; end

if ~exist('foldername','var')
    foldername = cd(cd('./'));
end
    files = list_h5_files(foldername);

for ii = 1:length(files)
    [folder,fname] = fileparts(files(ii).path);
    tmp = list_files(folder,sprintf('%s*eif.mat',fname));
    if isempty(tmp) | FORCE
        tmp = process_eif_file(files(ii).path,[],C);
    else
        fprintf(1,'Skipping file %s\n',files(ii).path)
    end
end

%% Merge all files
files = list_files(foldername,'*_eif.mat','all',{'trash'});
if ~isempty(files)
for ii = 1:length(files)
     data(ii) = load(files{ii});
end

tmp = [data.capacitance];
var.C = [tmp.C];
var.Vrest = [tmp.Vrest];
tmp = [data.eif];
var.param = vertcat(tmp.param);
tmp = [data.reif];
%%
ref_t = [];
ref_param = [];
for ii = 1:length(tmp)
    ref_t = [ref_t, tmp(ii).windows(1:end-1)];
%     ref_t = [ref_t, tmp(ii).windows(1:end-1)+diff(tmp(ii).windows)/2];
    ref_param = vertcat(ref_param,tmp(ii).param);
end
%%
var.ref_param = ref_param;
var.ref_t = ref_t;

if plotvar
    
    figure('visible','on');clf
    ax(1) = subplot(2,5,1); % Capacitance
    [count,bins] = hist(var.C);
    tmp = bar(bins,count);
    set(tmp,'facecolor','k','facecolor','k')
    tmp = xlim;
    xlim([min(tmp)-diff(tmp)/2,max(tmp)+diff(tmp)])
    xlabel('Capacitance (pF)')
    ylabel('Counts')
    
    label = {'Time constant (\tau_m)','Resting potential (E_m)', 'Spike width (/delta_T)','Spike threshold (V_t)'};
    for ii = 1:size(var.param,2) % eif parameters
        ax(ii+1) = subplot(2,5,ii+1);
        
        [count,bins] = hist(var.param(:,ii));
        tmp = bar(bins,count);
        set(tmp,'facecolor','k','facecolor','k')
        tmp = xlim;
        xlim([min(tmp)-diff(tmp)/2,max(tmp)+diff(tmp)])
        xlabel(label{ii})
        ylabel('Counts')
    
    end
    ax(2+size(var.param)) = subplot(2,5,ii+2); 
    [count,bins] = hist(var.Vrest);
    tmp = bar(bins,count);
    set(tmp,'facecolor','k','facecolor','k')
    tmp = xlim;
    xlim([min(tmp)-diff(tmp)/2,max(tmp)+diff(tmp)])
    xlabel('Estimated V_{rest} (mV)')
    ylabel('Counts')
    for ii = 1:size(var.param,2) 
        ax(ii+2+size(var.param,2)) = subplot(2,5,ii+2+size(var.param,2)); 
        plot(var.ref_t,var.ref_param(:,ii),'ko')
        hold on
        [~,edges] = max(arrayfun(@(x)length(x.windows),[data.reif]));
        edges = [data(edges).reif.windows,max(data(edges).reif.windows)+50]
        [x,y,ys] = binSamples(var.ref_t,var.ref_param(:,ii),edges);
        errorbar(x,y,ys)
        
        ylabel(label{ii})
        xlabel('time from spike (ms)')
    
    end
end
else
    var = [];
end
keyboard