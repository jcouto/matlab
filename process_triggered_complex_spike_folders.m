function process_triggered_complex_spike_folders(folders, ...
    spikelet_number,CSTHRESH,CSWINDOW,...
    MAX_SPIKELET,MAX_WAVES_PER_TRIAL,fromScratch)
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
%     - MAX_WAVES_PER_TRIAL consider only a max number of waves per trial.
%     - fromScratch (1/0) reprocess all raw data files.

figuresVisible='on';
CSLENGTH = true; % Plot the length of the complex spikes

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

if ~exist('MAX_WAVES_PER_TRIAL','var')
    MAX_WAVES_PER_TRIAL = 0; % Number of waveforms to use from
    % each trial. (0 = means use all)
end

if ~exist('fromScratch','var')
    fromScratch = 1; % Do from scratch
end
% For the captions
N_WAVES_parentesis = '';
if MAX_WAVES_PER_TRIAL
    N_WAVES_parentesis = sprintf([' (Including only a maximum of %d CS ' ...
        'waves per trial.)'],MAX_WAVES_PER_TRIAL);
end
RANDOM = false;
cwd = pwd;
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    expName = {'unknown'};
end
expName = expName{end};
% Merge complex spikes per folder and save them that way.
% if ~exist('spikelet_number','var')
all_cs = [];
for ii = 1:length(folders)
    foldername = cd(cd(folders{ii}));
    cd(folders{ii})
    files = list_files('.','*_cs.mat');
    h5files = list_h5_files('.');
    clear cs
    count = 1;
    for jj = 1:length(h5files)
        A = cellfun(@(x)~isempty(x),strfind(files, ...
            h5files(jj).path(1:end-3)));
        clear tmp
        if (sum(A) && ~fromScratch)
            tmp = load(files{A});
        else
            tmp = extract_triggered_complex_spikes(h5files(jj).path);
        end
        if ~isempty(tmp) & ~isempty(tmp.cs_waves)
            cs(count) = tmp;
            count = count+1;
            N = size(cs(end).cs_waves,1);
            if (length(unique(cs(end).stim_freq))>20)
                MAX_WAVES_PER_TRIAL = 0;
                disp('Using all waveforms. It is Noisy stim.')
                RANDOM = true;
            else
                disp(unique(cs(end).stim_freq))
            end
            
            if MAX_WAVES_PER_TRIAL
                if MAX_WAVES_PER_TRIAL < size(cs(end).cs_waves,1);
                    N = MAX_WAVES_PER_TRIAL;
                end
            end
            cs(end).cs_waves = cs(end).cs_waves(1:N,:);
            cs(end).stim_freq = cs(end).stim_freq(1:N);
            cs(end).cs = cs(end).cs(1:N);
            cs(end).iCSi = cs(end).iCSi(1:N);
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
                tmp = cellfun(@(x)x(unique_n(kk)),...
                    trial_isi_spikelets{jj}(trial_n_spikelets{jj}==unique_n(kk)));
                
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
    if CSLENGTH
        cslength = cellfun(@(x)x(N),cs_spikelets(n_spikelets==N));
        edges = min(cslength):0.3:max(cslength);
        bins = histc(cslength,edges);
        bar(edges,bins,'facecolor','k')
        plot([1,1].*mean(cslength),ylim,'--','color',cc(find(unique_n==N),:),'linewidth',2)
        xlabel(sprintf('Complex spike length (spikelet # %d) ',mode(n_spikelets)))
    else
        interspikelet = cellfun(@(x)x(N),isi_spikelets(n_spikelets==N));
        edges = min(interspikelet):0.3:max(interspikelet);
        bins = histc(interspikelet,edges);
        bar(edges,bins,'facecolor','k')
        plot([1,1].*mean(interspikelet),ylim,'--','color',cc(find(unique_n==N),:),'linewidth',2)
        xlabel(sprintf('ISI jitter (spikelet # %d) ',mode(n_spikelets)))
    end
    
    % Find experiment and folder name
    
    tmp = regexp(foldername,expName,'split');
    
    tmp = tmp{end};
    tmp(tmp =='/') = '_';
    appendix = sprintf('%s%s',expName,tmp);
    cs_diff = cellfun(@(x)diff(x),cs_timestamps,'uniformoutput',0);
    figName = sprintf('%s_cs_all.pdf',appendix);
    dataName = sprintf('%s_cs_all.mat',appendix);
    appendix(appendix =='_') = '-';
    if CSLENGTH
        Cstr = sprintf(['C - Histogram of the length of the Complex Spike measured ', ...
            'from the time of the last spikelet for CS with %d spikelets.'],N);
    else
        Cstr = sprintf(['C - Histogram of ',...
            'the jitter in the last spikelet for waveforms with %d spikelets. '],N);
    end
    
    caption = sprintf(['Experiment %s. Analysis of %d waveforms%s. Mean intercomplex spike interval: ',...
        '%3.2f s. A - Interspikelet interval versus the timing of the ',...
        'spikelet from the complex spike onset. Color codes for the max',...
        ' number of spikelets in the CS waveform. B - Histogram of the ',...
        'number of spikelets in each complex spike. %s',...
        'D - Complex spike raw voltage traces. Color code is the same as in A.',...
        ' E - For mean stimulus frequencies above 1Hz and each waveform with a',...
        ' particular number of spikelets and each trial; plot of the latency of',...
        ' the last spikelet versus the time from the beginning of the trial.'],...
        appendix,size(out.cs_waves,1),N_WAVES_parentesis,mean(cs_iCSi),Cstr);
    set(fig,'paperposition',[0,0,18,10],'papersize',[18,10])
    print(fig,'-dpdf',figName)
    printFigWithCaption(figName,caption,1)
    %%
    %Save data
    save(dataName,'-struct','out')
    cd(cwd)
    all_cs = [all_cs,cs];
end %% Done processing data from all folders
close all
% Plot a figure with the evolution of the experiment. I.e: "Actual"
% time versus the complex spike latency and stimulation frequency.
if ~MAX_WAVES_PER_TRIAL
    cc = setFigureDefaults();
    fig = figure('name','Experiment evolution','visible',figuresVisible);
    clear ax ann
    cs = all_cs;
    clear all_cs
    ax(1) = axes('position',[0.1,0.15,0.75,0.7]);
    ylabel('Complex spike duration')
    xlabel('Time (s)')
    % ax(2) = axes('position',[0.1,0.6,0.8,0.3]);
    % ylabel(sprintf('# detected spike(lets) in first %d ms',CSWINDOW))
    % xlabel('Time (s)')
    [~,recidx] = sort([cs.rec_date]);
    
    cs_abstime = [];
    cs_nspk = [];
    cs_length = [];
    stim_freq = [];
    % Reanalizing the data....not very efficient...
    for ii = recidx
        tmpstart = cs(ii).rec_date*24*60*60;
        
        cs_abstime = [cs_abstime;sort(cs(ii).cs + tmpstart)];
        % in seconds
        
        for jj = 1:size(cs(ii).cs_waves,1)
            tmpspks = cs(ii).twave(argfindpeaks(cs(ii).cs_waves(jj,:),CSTHRESH,...
                int32(0.6./diff(twave(1:2)))));
            cs_nspk = [cs_nspk, length(tmpspks<CSWINDOW)];
            cs_length = [cs_length, tmpspks(find(tmpspks<CSWINDOW,1,'last'))];
        end
        % is this this trial a unique frequency?
        tmp = unique(cs(ii).stim_freq);
        if length(tmp) == 1
            % freq, start,stop
            stim_freq = [stim_freq; tmp, tmpstart,cs(ii).rec_dur];
        end
    end
    cs_abstime = cs_abstime - cs(recidx(1)).rec_date*24*60*60;
    counter=1;
    for ii = unique(cs_nspk)
        p = plot(cs_abstime(cs_nspk == ii),cs_length(cs_nspk == ii),'ko', ...
            'markerfacecolor',cc(counter,:),'markersize',2);
        counter = counter+1;
    end
    axis tight
    if ~isempty(stim_freq)
        stim_freq(:,2) = stim_freq(:,2) - cs(recidx(1)).rec_date* ...
            24*60*60;
        for ii = 1:size(stim_freq,1)
            plot(stim_freq(ii,2)*[1,1],ylim,'r--')
            plot((stim_freq(ii,2)+stim_freq(ii,3))*[1,1],ylim,'k--')
            text(stim_freq(ii,2),max(ylim)*(0.8+randn(1)/5),...
                sprintf('%2.1fHz',stim_freq(ii,1)), ...
                'verticalalignment','top','horizontalalignment', ...
                'left')
        end
    end
    files = list_files('.','*_exp_evolution*.pdf',1);
    figName = sprintf('%s_exp_evolution_%02d.pdf',expName,length(files));
    
    title(sprintf(['Experiment %s - Complex spikes waveform length ']))
    set(fig,'paperposition',[0,0,18,8],'papersize',[18,7])
    print(fig,'-dpdf',figName)
end


%% Analyse the jitter of a particular spikelet.
if exist('spikelet_number','var')
    files = {};
    for ii = 1:length(folders)
        try
            tmp = list_files(folders{ii},'*_cs_all.mat');
        catch
            keyboard
        end
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
    if CSLENGTH
        for pp = 1:length(spikelets_time)
            try
                cslength(pp) = spikelets_time{pp}(end);
            catch
                cslength(pp) = nan;
            end
        end
    end
    
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
    
    files = list_files('.','*_cs_summary*.pdf',1);
    figName = sprintf('%s_cs_summary%02d.pdf',expName,length(files));
    if MAX_WAVES_PER_TRIAL
        figName = sprintf('%s_max%02dwaves_cs_summary%02d.pdf', ...
            expName, MAX_WAVES_PER_TRIAL,...
            length(files));
    end
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
    
    iCSi_BIN = 0.15;
    edges = (0.1:iCSi_BIN:10);
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
    if CSLENGTH
        plot(iCSi(analysis_idx),...
            cellfun(@(x)x(spikelet_number),...
            spikelets_time(analysis_idx)),...
            'ko','markeredgecolor','none','markerfacecolor',[.5,.5,.5],...
            'markersize',3)
    else
        plot(iCSi(analysis_idx),...
            cellfun(@(x)x(spikelet_number),...
            spikelets_isi(analysis_idx)),...
            'ko','markeredgecolor','none','markerfacecolor',[.5,.5,.5],...
            'markersize',3)
    end
    if length(intervals)>8
        temp = randperm(length(intervals));
        interv_idx = sort(temp(1:5));
    else
        interv_idx = 1:length(intervals);
    end
    iicount = 1;
    for ii = interv_idx
        idx = (iCSi>=(intervals(ii)) &...
            iCSi<(intervals(ii)+iCSi_BIN))' &...
            analysis_idx;
        axes(ax(1))
        try
            plot(twave,cs_waves(idx,:)+offset*(iicount-1),'color',cc(iicount,:))
        catch
            disp('vectors not the same length?')
            keyboard
        end
        axes(ax(2))
        plot(cellfun(@(x)x(spikelet_number),spikelets_time(idx)),...
            cellfun(@(x)x(spikelet_number),spikelets_isi(idx)),...
            'ko','markeredgecolor',cc(iicount,:),'markerfacecolor','none')
        axes(ax(4))
        if CSLENGTH
            plot(iCSi(idx),cellfun(@(x)x(spikelet_number),spikelets_time(idx)),...
                'ko','markeredgecolor',cc(iicount,:),'markerfacecolor','none')
        else
            plot(iCSi(idx),cellfun(@(x)x(spikelet_number),spikelets_isi(idx)),...
                'ko','markeredgecolor',cc(iicount,:),'markerfacecolor','none')
        end
        fid.cs_waves = [fid.cs_waves;cs_waves(idx,:)];
        fid.interCSinterval = [fid.interCSinterval;iCSi(idx)];
        
        iicount = iicount+1;
    end
    fid.spikeletNumber = spikelet_number;
    fid.all_iCSi = iCSi(analysis_idx);
    fid.spikelet_isi = cellfun(@(x)x(spikelet_number),spikelets_isi(analysis_idx))';
    fid.cs_length = cellfun(@(x)x(spikelet_number),spikelets_time(analysis_idx))';
    fid.expName = expName;
    fid.folders = folders;
    fid.random_stim=RANDOM;
    
    %%%%%%%%%%%%%%%%%%%%%
    axes(ax(4))
    
    G = {};
    F = {};
    warning off
    s1 = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0, 0.001, 2],...
        'Upper',[10, 1, 10]);
    ff{1} = fittype('a * (exp( -x / b )) + c',...
        'options',s1);
    s2 = fitoptions('Method','NonlinearLeastSquares','Lower',[-1,-20],...
        'Upper',[1, 20]);
    ff{2} = fittype('b * x + c','options',s2);
    fitidx = {find(fid.all_iCSi>0),find(fid.all_iCSi<=5)};
    fitname = {'Exponential','Linear'};
    if CSLENGTH
        iCSi_BIN = 0.2;
        edges = (0.2:iCSi_BIN:15)-0.05;
        [mean_x,mean_y,std_y,nn] = binSamples(fid.all_iCSi,fid.cs_length,edges);
        tmp_maxxaxis = 5;
        tmp_maxyaxis = 14;
        xlabel(sprintf('Complex spike length (spikelet # %d) ',spikelet_number))
        tmpx = fid.all_iCSi;
        tmpy = fid.cs_length;
        for f = 1:length(ff)
            if length(fitidx{f}) < 10
                fitidx{f} = find(fid.all_iCSi>0);
            end
            disp('Not enough points for the linear fit.')
        end
    else
        iCSi_BIN = 0.2;
        edges = (0.2:iCSi_BIN:7)-0.05;
        [mean_x,mean_y,std_y,nn] = binSamples(fid.all_iCSi,fid.spikelet_isi,edges);
        tmp_maxxaxis = 7;
        tmp_maxyaxis = 14;
        xlabel(sprintf('Inter-spikelet interval (spikelet # %d) ',spikelet_number))
        tmpx = fid.all_iCSi;
        tmpy = fid.spikelet_isi;
    end
    for f = 1:length(ff)
%         try
            [F{f},G{f}] = fit(tmpx(fitidx{f}),tmpy(fitidx{f}),ff{f});
%         catch
%             disp('Couldnt fit like this..')
%             whos
%             keyboard
%         end
    end
    
    warning on
    fid.fits = F;
    fid.gofits = G;
    
    tmp = F{1}; exp_f_val = tmp(edges);
    exp_pred_error = predint(F{1},edges,0.95);
    tmp = F{2}; lin_f_val = tmp(edges);
    fit_str = '';
    for ff = 1:2
        tmp = G{ff}
        fit_str = sprintf('%s %s fit, rsquared is %2.3f.',...
            fit_str,fitname{ff},tmp.rsquare);
    end
    idx = ~isnan(mean_y);
    errorbar(mean_x(idx),mean_y(idx),std_y(idx),'k',...
        'markerfacecolor',[.5,.5,.5])
    axis tight
    plot(edges,exp_f_val,'k-','linewidth',1.3)
    plot(edges,lin_f_val,'--','linewidth',1.3,'color',cc(2,:))
    plot(edges,exp_pred_error,'--','color',cc(1,:))
    
    if ~isempty(mean_y(find(nn>5,1,'last')))
        xlim([0.2,mean_y(find(nn>5,1,'last'))])
    else
        xlim([0.2,tmp_maxxaxis])
    end
    if max(ylim)>tmp_maxyaxis
        ylim([min(ylim),tmp_maxyaxis])
    end
    %%%%%%%%%%%%%%%%%%%%%
    axes(ax(1))
    axis tight
    xlim([min(xlim),15])
    set(ax(1),'ycolor','w')
    fname = [];
    for ii = 1:length(folders);fname = [fname,',',folders{ii}];end
    if CSLENGTH
        Cstr = sprintf(['C - Length of the Complex Spike (measured ', ...
            'from the time of the last spikelet for CS with %d ', ...
            'spikelets) versus the time from the last CS.%s '],N,fit_str);
    else
        Cstr = sprintf(['C - Inter-spikelet interval',...
            ' versus the time from the last CS.%s '],fit_str)
    end
    
    caption = sprintf(['Experiment %s - considering %d waveforms%s. A - Interspikelet interval versus the timing of the ',...
        'spikelet from the complex spike onset for the %d spikelet. Color codes for the ',...
        'inter-complex spike interval. B - Histogram of the ',...
        'intercomplex spike intervals for CS with %d spikelets. %s',...
        'D - Complex spike raw voltage traces. Color code is the ',...
        'same as in A. Used folders: %s'],...
        expName,length(find(analysis_idx)), ...
        N_WAVES_parentesis, ...
        spikelet_number,spikelet_number,Cstr,fname);
    
    set(fig,'paperposition',[0,0,18,10],'papersize',[18,10])
    print(fig,'-dpdf',figName)
    printFigWithCaption(figName,caption,1)
end
