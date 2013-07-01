function [units,units_spikes,tspkwave,units_wavelets,cluster_tree,cluster_input,cluster_temp] = process_extracellular_spikes(filename,entity_id,electrode_mapping)
% Extracts the spkwaves from a file using Quiroga's algorithms.
% In short the function follows this procedure:
%   - Filters raw traces
%   - Extracts the spkwaves with threshold crossing and dead time.
%   - Extracts the wavelet coefficients from each spike waveform.
%   - Runs SuperParamagnetic clustering on the wavelet coefficients.
%   - Analyses the temperature of the cluster tree and assigns the spkwaves to clusters.
%   - Saves spkwaves to a matlab file.

if ~exist('detection_mode','var')
    detection_mode = 'both';
end

tdead = 1.5;
tpost = 1;
tpre = 1;

[ent,info]=open_h5_trace(filename);

t = linspace(0,info.tend,length(ent(entity_id(1)).data(:)));
spks = cell(length(entity_id),1);
spkwaves = cell(length(entity_id),1);
threshold = cell(length(entity_id),1);
fprintf(1,'Extracting spkwaves from entity %d.\n',0);
for k = 1:length(entity_id)
    fprintf(1,'\b\b\b%d.\n',entity_id(k));
    V = ent(entity_id(k)).data(:);
    Vf = filter_data(V);
    
    [spk,spkwave,tspkwave,spk_idx,thresh] = extract_extracellular_spikes(Vf, [], t,tpre,tpost,tdead,detection_mode);
    spkwaves{k} = [spkwaves{k};spkwave];
    threshold{k} = [threshold{k},thresh];
    spks{k} = [spks{k},spk];
end
clear('V','Vf','ent')

fprintf(1,'Doing spike sorting ');
units = cell(1,length(spks));
units_spikes = cell(1,length(spks));
units_wavelets =  cell(1,length(spks));
cluster_input = cell(1,length(spks));
cluster_tree = cell(1,length(spks));
min_spk = 20;
discarded_units = cell(1,length(spks));
discarded_spkwaves = cell(1,length(spks));
for ii = 1:length(spks)
    fprintf(1,'.');
    [classes,cluster_tree{ii},cluster_input{ii},inspk] = compute_spike_sorting_classes(spkwaves{ii});
    
    tmp_units = {};
    tmp_shapes = {};
    tmp_wavelets = {};
    % Assign spkwaves to classes
    for jj = 1:length(classes)-1 % Last class are discarded spkwaves...
        if length(classes{jj})>min_spk
            tmp_units = [tmp_units,spks{ii}(classes{jj})];
            tmp_shapes = [tmp_shapes,spkwaves{ii}(classes{jj},:)];
            tmp_wavelets = [tmp_wavelets,inspk(classes{jj},:)];
        end
    end
    if isempty(classes{end})
        tmp_wavelets{end+1} = zeros(0,size(inspk,2));
    else
        tmp_wavelets{end+1} = inspk(classes{end},:);
    end
    
    
    discarded_units{ii} = spks{ii}(classes{end});
    discarded_spkwaves{ii} = spkwaves{ii}(classes{end},:);
    units{ii} = tmp_units;
    units_spikes{ii} = tmp_shapes;
    units_wavelets{ii} = tmp_wavelets;
end
[cluster_temp] = find_tree_temperature(cluster_tree{ii},cluster_input{ii});
fprintf(1,'\n');


[tmppath,file_name] = fileparts(filename);
visibilityflag = 'on';
if ~exist('plot_shapes','var')
    visibilityflag = 'off';
end
cc = setFigureDefaults();
smooth_bin = 0.1; % in seconds
smooth_window = 30;
tedges = (0:smooth_bin:info.tend);
isiedges = 0.001*10.^((1:35+1)/9); % use a log scale
for ii = 1:length(units)
    fig = figure(1);set(fig,'Visible',visibilityflag);clf
    
    ax(1) = axes('position',[0.1,0.1,0.4,0.4]);hold on
    ax(2) = axes('position',[0.55,0.1,0.4,0.4]);hold on
    ax(3) = axes('position',[0.1,0.65,0.6,0.25]);hold on
    ax(4) = axes('position',[0.75,0.65,0.2,0.25]);hold on
    if ~isempty(discarded_spkwaves{ii})
        axes(ax(1))
        plot(tspkwave,discarded_spkwaves{ii},'color',[0.3,0.3,0.3],'linewidth',0.5)
        axes(ax(4))
        plot(units_wavelets{ii}{end}(:,1),units_wavelets{ii}{end}(:,2),'ko',...
            'markerfacecolor','none','markeredgecolor',[0.3,0.3,0.3],'markersize',3)
    end
    for kk = 1:length(units{ii})
        axes(ax(1))
        plot(tspkwave,units_spikes{ii}{kk},'color',cc(kk,:),'linewidth',0.5)
        axes(ax(3))
        binned_sp = histc(units{ii}{kk},tedges);
        w = gausswin(smooth_window*smooth_bin);
        y = filter(w,1,binned_sp./smooth_bin);
        plot(tedges,y,'color',cc(kk,:))
        axes(ax(2))
        binned_isi = histc(diff(units{ii}{kk}),isiedges);
        h = bar(isiedges,binned_isi,'histc');
        set(h,'edgecolor',cc(kk,:),'facecolor','none')
        
    end
    axes(ax(2))
    axis tight
    set(ax(2),'xscale','log','xtick',[0.003,0.01,0.1,1,10],'xticklabel',...
        [3,10,100,1e3,10e3],'xlim',[0,10])
    
    
    for kk = 1:length(units{ii})
        axes(ax(1))
        plot(tspkwave,mean(units_spikes{ii}{kk}),'color',cc(kk,:)*0.6,'linewidth',2)
        axes(ax(4))
        plot(units_wavelets{ii}{kk}(:,1),units_wavelets{ii}{kk}(:,2),'ko',...
            'markerfacecolor','none','markeredgecolor',cc(kk,:)*0.6,'markersize',3)
    end
    axes(ax(1)),axis tight
    axes(ax(4)),axis tight
    %[temp] = find_tree_temperature(cluster_tree{ii},cluster_input{ii},[],15,[],[],1);
    %set(ax(4),'xscale','log','yscale','log')
    set(fig,'paperposition',[0,0,18,15],'papersize',[18,15],'paperunits','centimeters')
    print(fig,'-dpdf',sprintf('%s/%s_sorted_entity%3.0d.pdf',tmppath,file_name,entity_id(ii)))
end

save(sprintf('%s/%s_sorted_spikes.mat',tmppath,file_name),'filename','units','units_spikes',...
    'tspkwave','units_wavelets','cluster_tree','cluster_input','cluster_temp');

