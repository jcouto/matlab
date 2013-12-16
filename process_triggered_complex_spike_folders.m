function process_triggered_complex_spike_folders(folders, spikelet_number,CSTHRESH,CSWINDOW,MAX_SPIKELET,fromScratch)
% Extracts the complex spikes from the given folders and produces a
% figure and saves the traces.
% Input parameters:
%     - folders is a cell array with the names of the folders
%     - spikelet number is the number of the spikelets to be
%     accounted for.
%     - CSTHRESH is the threshold for spikelet detection.
%     - CSWINDOW is the window for spikelet detection
%     - MAX_SPIKELET is another max that can be used to further
%     sort the spikes. It represents the maximum of the spikelet
%     that is being analysed.
%     - fromScratch (1/0) reprocess all raw data files.

figuresVisible='on';

if ~exist('folders','var')
    folders = dir();
    folders = {folders([folders.isdir]).name};
    folders = folders(3:end);
end
if ~exist('CSTHRESH','var')
    CSTHRESH = -40; % threshold for spikelet detection
end
if ~exist('CSWINDOW','var')
    CSWINDOW = 15; % window for spikelet detection
end
if ~exist('MAX_SPIKELET','var')
    MAX_SPIKELET = 50; % max time for principal spikelet
end
if ~exist('fromScratch','var')
    fromScratch = 1; % Do from scratch
end
cwd = pwd;
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    expName = {'unknown'};
end
expName = expName{end};
% Merge complex spikes per folder and save them that way.
% if ~exist('spikelet_number','var')
for ii = 1:length(folders)
    foldername = cd(cd(folders{ii}));
    cd(folders{ii})
    files = list_files('.','*_cs.mat');
    h5files = list_h5_files('.');
    clear cs
    count = 1;
    for jj = 1:length(h5files)
        A = cellfun(@(x)~isempty(x),strfind(files,h5files(jj).path(1:end-3)));
        if (sum(A) && ~fromScratch)
            tmp = load(files{A});
        else
            tmp = extract_triggered_complex_spikes(h5files(jj).path);
        end
        if ~isempty(tmp.cs_waves)
            cs(count) = tmp;
            count = count+1;
        end
    end
    
    trial_cs_waves = {cs.cs_waves};
    cs_waves = vertcat(cs.cs_waves);
    twave = cs(1).twave;
    nstim_bursts = vertcat(cs.nstim_bursts);
    stim_freq = horzcat(cs.stim_freq);
    cs_timestamps = {cs.cs};
    cs_iCSi = vertcat(cs.iCSi);

    % Compute the cs_spikelets on a trial basis
    trial_cs_spikelets = cell(1,length(trial_cs_waves));
    trial_isi_spikelets = cell(1,length(trial_cs_waves));
    trial_n_spikelets = cell(1,length(trial_cs_waves));
    
    for jj = 1:size(trial_cs_waves,2)
        trial_cs_spikelets{jj} = cell(1,size(trial_cs_waves{jj},1));
        for kk = 1 :size(trial_cs_waves{jj},1)
            idx = argfindpeaks(trial_cs_waves{jj}(kk,:),CSTHRESH,...
                               int32(0.6./diff(twave(1:2))));
            cs_lim = find(twave>0 & twave<CSWINDOW);
            tmp = twave(idx(idx>min(cs_lim) & idx<max(cs_lim)));
            trial_cs_spikelets{jj}{kk} = tmp;
        end
        trial_isi_spikelets{jj} = cellfun(@(x)diff([0,x]),trial_cs_spikelets{jj},...
                                          'uniformoutput',0);
        trial_n_spikelets{jj} = cellfun(@length,trial_cs_spikelets{jj});
    end
    % Extract the spikelets timing
    %    CSWINDOW = 10; % window to search for spikelets!
    cs_spikelets = cell(1,size(cs_waves,1));
    for kk = 1 :size(cs_waves,1)
        idx = argfindpeaks(cs_waves(kk,:),CSTHRESH,int32(0.6./diff(twave(1:2))));
        cs_lim = find(twave>0 & twave<CSWINDOW);
        tmp = twave(idx(idx>min(cs_lim) & idx<max(cs_lim)));
        cs_spikelets{kk} = tmp;
    end
    % % Interspikelet interval versus spikelet time
    isi_spikelets = cellfun(@(x)diff([0,x]),cs_spikelets,'uniformoutput',0);
    max_n = max(cellfun(@length,cs_spikelets));
    n_spikelets = cellfun(@length,cs_spikelets);
    unique_n = unique(n_spikelets);
    unique_n(unique_n==0) = [];
    for jj = 1:length(unique_n)
        idx = find(n_spikelets==unique_n(jj));
        X{jj}=cell(1,max_n);
        Y{jj}=cell(1,max_n);
        for kk = 1:max_n
            for pp = idx
                try
                    X{jj}{kk} = [X{jj}{kk},cs_spikelets{pp}(kk)];
                    Y{jj}{kk} = [Y{jj}{kk},isi_spikelets{pp}(kk)];
                end
            end
        end
    end
    out.cs_waves = cs_waves;
    out.twave = twave;
    out.stim_bursts = nstim_bursts;
    out.stim_freq = stim_freq;
    out.cs_timestamps = cs_timestamps;
    out.cs_iCSi = cs_iCSi;
    out.spikelets_time = cs_spikelets;
    out.spikelets_isi = isi_spikelets;
    out.spikelets_count = n_spikelets;
    out.isolated_spikelets_time = X;
    out.isolated_spikelets_isi = Y;

    %% Figure for the folder!
    cc = setFigureDefaults();
    fig = figure('name',['Folder ',folders{ii}],'visible',figuresVisible);
    
    ax(1) = axes('position',[0.1,0.35,0.8,0.25]);
    ylabel('Voltage (mV)')
    xlabel('Time (ms)')
    ax(2) = axes('position',[0.1,0.7,0.2,0.2]);
    ylabel('Inter-spikelet (ms)')
    xlabel('Spikelet time (ms)')
    ax(3) = axes('position',[0.4,0.7,0.2,0.2]);
    ylabel('# of Complex Spikes')
    xlabel('# Spikelets')
    ax(4) = axes('position',[0.7,0.7,0.2,0.2]);
    ylabel('Count')
    xlabel(sprintf('ISI jitter (spikelet # %d) ',mode(n_spikelets)))
    ax(5) = axes('position',[0.1,0.1,0.8,0.18]);
    ylabel({'Last isi (ms)'})
    xlabel({'Time (s) from trial start ','(only for f>1Hz)'})
    set(ax(2:end-1),'ticklength',[0.025,0.03])
    label = 'DABCE';
    for kk = 1:length(ax)
        p = get(ax(kk),'position');
        ann(kk) = annotation('textbox',p+[-0.06,0.06,0,0],'string',label(kk));
    end
    set(ann,'verticalalignment','top',...
        'horizontalalignment','left',...
        'color','k','edgecolor','none',...
        'fontsize',10,'fontweight','bold')
    
    % Interspikelet interval versus spikelet time
    for kk = 1:length(unique_n)
        axes(ax(2))
        for jj = 1:max_n
            plot(X{kk}{jj},Y{kk}{jj},'ko',...
                'markeredgecolor',cc(kk,:),...
                'markerfacecolor','none','markersize',2.5)
        end
        axes(ax(1))
        % Complex spike shape
        plot(twave,cs_waves(n_spikelets==unique_n(kk),:),'color', ...
             cc(kk,:),'linewidth',0.6)
        axes(ax(5))
        text_location = [];
        for jj = 1:length(trial_isi_spikelets)
            if ~isempty(find(trial_n_spikelets{jj}==unique_n(kk)))
                try % Not needed...
                    tmp = cellfun(@(x)x(unique_n(kk)),...
                                  trial_isi_spikelets{jj}(trial_n_spikelets{jj}==unique_n(kk)));
                catch
                    keyboard
                end
                if ~isempty(tmp)
                    if ((1./mean(cs(jj).iCSi))>1) && length(tmp)>15
                        tmpx = out.cs_timestamps{jj}(trial_n_spikelets{jj}==unique_n(kk));
                        plot(tmpx,tmp,...
                             'k-o','color',cc(kk,:),'linewidth',0.6,'markersize',2,...
                             'markerfacecolor','none','markeredgecolor', ...
                             cc(kk,:))
                        if sum(((tmp(end)+1)<text_location) & ((tmp(end)-1)>text_location))
                            text_location = [text_location,tmp(end)+1];
                        else
                            text_location = [text_location,tmp(end)];
                        end
                        
                        text(tmpx(end),text_location(end),sprintf(['trial %d, %3.2fHz'], ...
                                                                  jj,1./mean(cs(jj).iCSi)),...
                             'verticalalignment','bottom',...
                             'horizontalalignment','left',...
                             'color',cc(kk,:))
                    end
                end
            end
        end
    end
    axes(ax(5))
    axis tight
    if max(ylim)>CSWINDOW
        ylim([min(ylim),CSWINDOW])
    end
    
    % Histogram of the number of spikelets
    axes(ax(3))
    bins = histc(n_spikelets,0:max_n);
    
    bar(0:max_n,bins,'facecolor','k')
    axis tight
    xlim([0,max(xlim)+1])
    % Jitter of the last spikelet
    axes(ax(4))
    N=mode(n_spikelets);
    interspikelet = cellfun(@(x)x(N),isi_spikelets(n_spikelets==N));
    edges = min(interspikelet):0.3:max(interspikelet);
    bins = histc(interspikelet,edges);
    bar(edges,bins,'facecolor','k')
    plot([1,1].*mean(interspikelet),ylim,'--','color',cc(find(unique_n==N),:),'linewidth',2)
    
    % Find experiment and folder name
    
    tmp = regexp(foldername,expName,'split');
    
    tmp = tmp{end};
    tmp(tmp =='/') = '_';
    appendix = sprintf('%s%s',expName,tmp);
    cs_diff = cellfun(@(x)diff(x),cs_timestamps,'uniformoutput',0);
    figName = sprintf('%s_cs_all.pdf',appendix);
    dataName = sprintf('%s_cs_all.mat',appendix);
    appendix(appendix =='_') = '-';
    caption = sprintf(['Experiment %s. Analysis of %d waveforms. Mean intercomplex spike interval: ',...
                       '%3.2f s. A - Interspikelet interval versus the timing of the ',...
                       'spikelet from the complex spike onset. Color codes for the max',...
                       ' number of spikelets in the CS waveform. B - Histogram of the ',...
                       'number of spikelets in each complex spike. C - Histogram of ',...
                       'the jitter in the last spikelet for waveforms with %d spikelets. ',...
                       'D - Complex spike raw voltage traces. Color code is the same as in A.',...
                       ' E - For mean stimulus frequencies above 1Hz and each waveform with a',...
                       ' particular number of spikelets and each trial; plot of the latency of',...
                       ' the last spikelet versus the time from the beginning of the trial.'],...
        appendix,size(out.cs_waves,1),mean(cs_iCSi),N);
    set(fig,'paperposition',[0,0,18,10],'papersize',[18,10])
    print(fig,'-dpdf',figName)
    printFigWithCaption(figName,caption,1)
    %%
    %Save data
    save(dataName,'-struct','out')
    cd(cwd)
end

%% Analyse the jitter of a particular spikelet.
if exist('spikelet_number','var')
    files = {};
    for ii = 1:length(folders)
        tmp = list_files(folders{ii},'*_cs_all.mat');
        files = horzcat(files,tmp);
        tmp = load(files{end});
        data(ii) = tmp;
        try;close(['Folder ',folders{ii}]);end
    end
    %% Plot the waveforms only for a particular number of spikelets and inter CS interval
    iCSi = vertcat(data.cs_iCSi);
    cs_waves = vertcat(data.cs_waves);
    twave = data(1).twave;
    spikelet_count = horzcat(data.spikelets_count);
    spikelets_time = horzcat(data.spikelets_time);
    spikelets_isi = horzcat(data.spikelets_isi);
    %%
    cc = setFigureDefaults();
    fig = figure('name','Summary CS','visible',figuresVisible);
    clear ax ann
    ax(1) = axes('position',[0.1,0.1,0.8,0.5]);
    ylabel('Voltage (mV)')
    xlabel('Time (ms)')
    ax(2) = axes('position',[0.1,0.7,0.2,0.2]);
    ylabel('Inter-spikelet (ms)')
    xlabel('Spikelet time (ms)')
    ax(3) = axes('position',[0.4,0.7,0.2,0.2]);
    ylabel('Counts (# CS)')
    xlabel(sprintf('Time from previous CS (s)'))
    
    ax(4) = axes('position',[0.7,0.7,0.2,0.2]);
    ylabel('inter-spikelet interval (ms)')
    xlabel(sprintf('Time from previous CS (s)'))
    set(ax(2:end),'ticklength',[0.025,0.03])
    label = 'DABC';
    for kk = 1:length(ax)
        p = get(ax(kk),'position');
        ann(kk) = annotation('textbox',p+[-0.06,0.06,0,0],'string',label(kk));
    end
    set(ann,'verticalalignment','top',...
        'horizontalalignment','left',...
        'color','k','edgecolor','none',...
        'fontsize',10,'fontweight','bold')
    
    %     MAX_SPIKELET = 6;
    %    MAX_SPIKELET = 50;
    files = list_files('.','*_cs_summary*.pdf',1);
    figName = sprintf('%s_cs_summary%02d.pdf',expName,length(files));
    matName = sprintf('%s_cs_summary%02d.mat',expName,length(files));

    fid = matfile(matName);
    fid.cs_waves = {};
    fid.twaves = twave;
    fid.interCSinterval = {};

    
    tmp = nan(size(spikelets_time));
    for ii = 1:length(spikelets_time)
        try
        tmp(ii) = spikelets_time{ii}(spikelet_number);
        end
    end
    analysis_idx = (spikelet_count==spikelet_number)&(tmp<MAX_SPIKELET); 

    iCSi_BIN = 0.2;
    edges = (0.15:iCSi_BIN:7)-0.05;
    axes(ax(3))

    bins = histc(iCSi(analysis_idx),edges);
    bar(edges,bins,'facecolor','k')
    intervals = edges(bins>3);
    
    offset = 50;
    try
        xlim([min(intervals),max(intervals)])
    catch
        disp('Cant set ylim... intervals empty or the same?')
    end
    axes(ax(4))
    plot(iCSi(analysis_idx),...
         cellfun(@(x)x(spikelet_number),...
                 spikelets_isi(analysis_idx)),...
         'ko','markeredgecolor','none','markerfacecolor',[.5,.5,.5],...
         'markersize',3)
    if length(intervals)>8
        temp = randperm(length(intervals));
        interv_idx = sort(temp(1:5));
    else
        interv_idx = 1:length(intervals);
    end
    iicount = 1;
    for ii = interv_idx
        idx = (iCSi>(intervals(ii)) &...
            iCSi<=(intervals(ii)+iCSi_BIN))' &...
             analysis_idx;
        axes(ax(1))
        plot(twave,cs_waves(idx,:)+offset*(iicount-1),'color',cc(iicount,:))
        axes(ax(2))
        plot(cellfun(@(x)x(spikelet_number),spikelets_time(idx)),...
             cellfun(@(x)x(spikelet_number),spikelets_isi(idx)),...
             'ko','markeredgecolor',cc(iicount,:),'markerfacecolor','none')
        axes(ax(4))
        plot(iCSi(idx),cellfun(@(x)x(spikelet_number),spikelets_isi(idx)),...
            'ko','markeredgecolor',cc(iicount,:),'markerfacecolor','none')
             fid.cs_waves = [fid.cs_waves;cs_waves(idx,:)];
        fid.interCSinterval = [fid.interCSinterval;iCSi(idx)];
        iicount = iicount+1;
    end
    iCSi_BIN = 0.2;
    edges = (0.2:iCSi_BIN:7)-0.05;
    analysis_idx = (spikelet_count==spikelet_number) & (tmp<MAX_SPIKELET);
    [mean_iCSi,mean_spikelet_isi,std_spikelet_isi] = binSamples(iCSi(analysis_idx),...
                                                      cellfun(@(x)x(spikelet_number),...
                                                      spikelets_isi(analysis_idx)),...
                                                      edges);

    axes(ax(4))
    idx = ~isnan(mean_spikelet_isi);
    errorbar(mean_iCSi(idx),mean_spikelet_isi(idx),std_spikelet_isi(idx),'k','markerfacecolor',[.5,.5,.5])
    xlim([min(mean_iCSi(idx)-0.2),max(mean_iCSi(idx))+0.2])
    axes(ax(1))
    axis tight
    xlim([min(xlim),15])
    set(ax(1),'ycolor','w')
    fname = [];for ii = 1:length(folders);fname = [fname,',',folders{ii}];end
    caption = sprintf(['Experiment %s - considering %d waveforms. A - Interspikelet interval versus the timing of the ',...
        'spikelet from the complex spike onset for the %d spikelet. Color codes for the ',...
        'inter-complex spike interval. B - Histogram of the ',...
        'intercomplex spike intervals for CS with %d spikelets. C - Inter-spikelet interval',...
        ' versus the time from the last CS. ',...
        'D - Complex spike raw voltage traces. Color code is the same as in A. Used folders: %s'],...
        expName,length(find(analysis_idx)),spikelet_number,spikelet_number,fname);
    
    set(fig,'paperposition',[0,0,18,10],'papersize',[18,10])
    print(fig,'-dpdf',figName)
    printFigWithCaption(figName,caption,1)
    %     axes(ax(1))
end
