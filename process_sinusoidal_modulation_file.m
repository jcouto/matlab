function [data,found_sinusoid] = process_sinusoidal_modulation_file(filename,nBins,threshold,detection_mode,holding_current,dendrite)
% This function processes a file that had a sinusoidal waveform injected.
% [data,found_sinusoid] = process_sinusoidal_modulation_file(filename,nBins,threshold,detection_mode,dendrite)
%
%

[files,kfiles] = list_h5_files;
[~,filename,ext] = fileparts(filename);
file = files(find(~cellfun(@isempty,(strfind({files.basename},filename)))));
found_sinusoid = 0;
if ~exist('nBins','var');nBins = 10;end
if ~exist('plotvar','var');plotvar = 1;end
if ~exist('threshold','var')||isempty(threshold);threshold = -25;end
if ~exist('holding_current','var');holding_current = [0,0];end
if ~exist('spike_alignment_mode','var');spike_alignment_mode = 'peak';end % peak or thresh
if ~exist('dendrite','var');dendrite = [];end
if ~exist('sortByTrial');sortByTrial=false;end

[tall,Vall,Iall,W,metadata,info,Vdall,Imall] = i_load_compensated_voltage(file,kfiles,holding_current,dendrite);
% Assumes the last thing is the sinusoids protocol (and that there is a
% tail just before the end).
tprot = unique(cumsum(metadata{1}(:,1)));
% tstart_sin = tprot(end-2);
% tend_sin = tprot(end-1);
% V0 = nanmean(Vall(tall<tprot(1)));

dt = tall(2)-tall(1);
% dVdt = (diff(Vall)*1.0e-3)./dt;

[spk,spkwave,tspkwave] = extract_spikes( Vall, threshold, tall, 1, 3, 4);
% Find Vm
tprot = unique(cumsum(metadata{1}(:,1)));
spont.Vm = nanmean(Vall(tall<tprot(1)));
spont.Vm_std = std(Vall(tall<tprot(1)));

% Find spontaneous firing
spont.spks = spk(spk <  tprot(1));
isi = diff(spont.spks);
spont.frate = 1./nanmean(isi);
spont.cv = std(isi)/nanmean(isi);
spont.spkwave = spkwave(spk < tprot(1),:);
spont.tspkwave = tspkwave;

% Spike parameters

[spont.spk_Vthr,spont.spk_height,spont.spk_slope,...
    spont.spk_width,spont.spk_ahdp] = extract_AP_features...
    (spont.spkwave,dt);

% Stimulation parameters

noiseidx = cellfun(@(x)find(x(:,2) == 2 | x(:,10) == 2),metadata,'uniformoutput',0);
sinidx = cellfun(@(x)find(x(:,2) == 3 | x(:,10) == 3),metadata,'uniformoutput',0);

stim.type = 'conductance';
if strcmp(W.units,'pA')
    stim.type = 'current';
end
tmpA = cellfun(@(x)~isempty(x),noiseidx);
tmpB = cellfun(@(x)~isempty(x),sinidx);
stim.mean = metadata{tmpA}(noiseidx{tmpA},3);  %stim nanmean
stim.sigma = metadata{tmpA}(noiseidx{tmpA},4); %noise standard deviation
stim.tau = metadata{tmpA}(noiseidx{tmpA},5); %noise time constant

stim.F = metadata{tmpB}(sinidx{tmpB},4) %sinusoidal modulation frequency
if stim.F > 2000
    keyboard
end
stim.m = metadata{tmpB}(sinidx{tmpB},3); %sinusoidal modulation amplitude
stim.tstart = tprot(end-2); %start of stim
stim.tend = tprot(end-1); %end of stim
stim.sinusoidal = @(t)stim.m*sin(2*pi*stim.F*t+0); %if it is a sinusoidal

% Estimate capacitance from noisy trace
if (stim.F > 0) && (abs(stim.m) > 0)
    found_sinusoid = 1;
end

spont.C = nan;
spont.Ce = [];
spont.var_I = [];

if (stim.sigma>0 || stim.tau>0)
    cidx = (tall>stim.tstart+2 & tall<stim.tend-2);
    [spont.C,spont.Ce,spont.var_I] = estimate_capacitance_from_noisy_trace...
        (tall(cidx),Vall(cidx),Iall(cidx),spont.Vm);
end
stim.Vm = mean(Vall(cidx));
stim.Vm_std = std(Vall(cidx));

int_factor = 3;

% Interpolates the spike waveforms to get the corrected spike time
stim_spk = find(spk>stim.tstart&spk<stim.tend);
[spk(stim_spk),spk_correction] = int_spikes(spk(stim_spk),spkwave(stim_spk,:),dt,int_factor,spike_alignment_mode);
% Values used to correct spike times
spont.spk_correction_mean = nanmean(spk_correction);
spont.spk_correction_std = nanmean(spk_correction);

% prepare to save data
data.spont = spont;
data.stim = stim;
data.spiketimes = spk;
data.spkwaves = spkwave;
data.tspkwaves = tspkwave;
data.threshold = threshold;
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
[~,filename] = fileparts(file.basename);
data.expName = expName{end};
data.filename = filename;
switch found_sinusoid
    case false
        isi = diff(spk(spk>stim.tstart&spk<stim.tend));
        data.cv = std(isi)/nanmean(isi);
        data.fr = 1./nanmean(isi);
        save(sprintf('no_modulation_data%s.mat',filename),'-struct','data');
    otherwise
    % Extract modulation
    [r0,r1,phi,conf,r,relspks,model] = compute_sinusoidal_modulation(...
        spk(spk>stim.tstart&spk<stim.tend)- stim.tstart, ...
        stim.F, nBins, stim.tend - stim.tstart, 1, 1);
    
    % Extract pseudo-trials from relative spikes.
    relspks_idx = [1,find(diff(relspks)<=0)+1];
    spk_mod = cell(length(relspks_idx)-1,1);
    for ii = 1:length(relspks_idx)-1
        spk_mod{ii}  = relspks(relspks_idx(ii):relspks_idx(ii+1)-1)*1.0e3; %in ms
    end
    %%% extract the sinusoidal modulation from the shuffled spike times
    shuffle_iter = 50;
    for ii=1:shuffle_iter
        % shuffle the spike times
        shuffled_spk = (spk(spk>stim.tstart&spk<stim.tend)- stim.tstart);
        isi = diff(shuffled_spk);
        idx = randperm(length(isi));
        shuffled_spk = cumsum(isi(idx));
        [r0_shuffled(ii),r1_shuffled(ii),phi_shuffled(ii)] = ...
            compute_sinusoidal_modulation(shuffled_spk, ...
            stim.F, nBins, stim.tend - stim.tstart, 1, 1);
        
    end
    
    data.r0_shuffled = nanmean(r0_shuffled);
    data.r1_shuffled = nanmean(r1_shuffled);
    data.phi_shuffled = nanmean(phi_shuffled);
    
    % Save data
    
    data.r0 = r0;
    data.r1 = r1;
    data.phi = phi;
    data.conf = conf;
    data.r = r;
    data.spk_mod = spk_mod;
    data.model = model;
    save(sprintf('sin_mod_data_%s.mat',filename),'-struct','data');
    %% Plot figure (for this file)
    fig = figure();
    cc = setFigureDefaults;
    ax(1) = axes('position',[0.1,0.8,.8,.15]);
    pnperiods = 8;
    
    tplot = spk(ceil(length(spk)/2)) - (pnperiods/2)./stim.F;%stim.tstart+(stim.tend-stim.tstart)/2;
    p = plot(tall(tall>tplot&tall<tplot+pnperiods./stim.F),...
        Vall(tall>tplot&tall<tplot+pnperiods./stim.F),'k'); % plot 5 periods
    set(p,'clipping','off')
    hold on
    
    % Plot spike times
    plot(repmat(spk(spk>tplot&spk<tplot+pnperiods./stim.F),2,1),...
        repmat(spont.Vm+[-3,2]',1,length(spk(spk>tplot&spk<tplot+pnperiods./stim.F))),...
        'color',cc(1,:),'linewidth',1.2)
    ylim(spont.Vm+[-10,+5])
    if dendrite
        plot(tall(tall>tplot&tall<tplot+pnperiods./stim.F),...
            Vdall(tall>tplot&tall<tplot+pnperiods./stim.F),'k','color',cc(3,:))
        ylim([min([spont.Vm-10,min(Vdall)]),max([spont.Vm+5,max(Vdall)])])
    end
    ylabel({'Membrane voltage (mV)','(spikes are trimmed)'})
    % Plot sinusoid
    ax(2) = axes('position',get(ax(1),'position'),...
        'yaxislocation','right','ycolor',cc(2,:));
    tstim = tall(tall>stim.tstart & tall<=stim.tend);
    plot(tall(tall>tplot&tall<tplot+pnperiods./stim.F),...
        stim.sinusoidal(tall(tall>tplot&tall<tplot+pnperiods./stim.F)-tstim(1)),...
        'color',cc(2,:),'linewidth',1.5) % plot 5 periods
    ylim([-3*stim.m,3*stim.m])
    xlabel('Time (s)')
    ylabel({'Sinusoidal modulation','current(pA)'})
    linkaxes(ax,'x')
    xlim([tplot,tplot+pnperiods./stim.F])
    % hold on
    % plot(tall(spkidx),Vall(spkidx),'ko')
    % ax1 = gca;
    % ax2 = axes();
    % set(ax2,'color','none','yaxislocation','right')
    % linkaxes([ax1,ax2],'x')
    % plot(tall,W.data,'go');
    % tstim = tall(tall>stim.tstart & tall<=stim.tend);
    % plot(tstim,stim.sinusoidal(tstim - tstim(1)),'r')
    ax(3) = axes('position',[0.1,0.45,.4,.3]);
    edges = linspace(0,1000./stim.F,length(r));
    x = linspace(0,1000./stim.F,nBins*10);
    hold on;
    plot([edges(1),edges(end)],r0+[0,0],'--','color',cc(2,:),'LineWidth',1.2);
    hndl = bar(edges, r, 0.8, 'histc');
    set(hndl,'FaceColor',[0.5,.5,.5],'linewidth',1);
    plot(x, r0+r1*sin(2*pi*stim.F*x/1000.0), '--','color',cc(3,:), 'LineWidth', 1.5);
    plot(x, r0+r1*sin(2*pi*stim.F*x/1000.0+phi),'color',cc(1,:), 'LineWidth', 1.2);
    axis([0,edges(end),min(r(1:end-1))-5,max(r)+5]);
    % Make a rasterplot triggered to the sinusoid.
    ax(4) = axes('position',[0.1,0.1,.4,.3],'ytick',[],'ycolor','w');
    [~,sortidx]=sort(cellfun(@(z)min(z),spk_mod));
    if sortByTrial
        plotRastergram(spk_mod(sortidx),'type','few');
    else
        plotRastergram(spk_mod,'type','few');
    end
    axis tight;xlim([0,1000./stim.F])
    xlabel('Time from sinusoidal cycle (ms)')
    linkaxes(ax([3,4]),'x')
    if ~isempty(spont.spks)
        ax(5) = axes('position',[0.6,0.45,.3,.25]);
        plot(spont.tspkwave(spont.tspkwave<1.5),spont.spkwave(:,(spont.tspkwave<1.5)),'k')
        plot(spont.tspkwave(spont.tspkwave<1.5),nanmean(spont.spkwave(:,(spont.tspkwave<1.5))),'color',cc(1,:))
        plot(xlim,nanmean(spont.spk_Vthr)*[1,1],'--','color',cc(2,:))
        text(max(xlim),nanmean(spont.spk_Vthr),sprintf('Threshold %3.1f mV',nanmean(spont.spk_Vthr)),...
            'VerticalAlignment','Bottom','HorizontalAlignment','Right','color',cc(2,:))
        plot(xlim,nanmean(spont.spk_height+spont.spk_Vthr)*[1,1],'--','color',cc(3,:))
        text(max(xlim),nanmean(spont.spk_height+spont.spk_Vthr),sprintf('Height %3.1f mV',nanmean(spont.spk_height)),...
            'VerticalAlignment','Bottom','HorizontalAlignment','Right','color',cc(3,:))
        xlabel('Time from spike peak (ms)')
        ylabel('V_m (mV)')
        axis tight
        ax(6) = axes('position',[0.6,0.1,.3,.25]);
        d1Xi = diff(spont.spkwave,1,2)./(dt*1.e3);
        plot(spont.spkwave(:,1:end-1)',d1Xi','k')
        plot(nanmean(spont.spkwave(:,1:end-1)),nanmean(d1Xi),'color',cc(1,:))
        axis tight
        plot(nanmean(spont.spk_Vthr)*[1,1],ylim,'--','color',cc(2,:))
        ylabel('dV/dt (mV.ms^{-1})')
        xlabel('V_m (mV)')
        ax(7) = axes('position',[0.75,0.2,.15,.15],'color','w');
        plot(spont.spkwave(:,1:end-1)',d1Xi','k')
        plot(nanmean(spont.spkwave(:,1:end-1)),nanmean(d1Xi),'color',cc(1,:))
        axis tight
        plot(nanmean(spont.spk_Vthr)*[1,1],ylim,'--','color',cc(2,:))
        ylabel('dV/dt (mV.ms^{-1})')
        xlabel('V_m (mV)')
        tmpmean=nanmean(d1Xi);
        a = find(nanmean(spont.spk_Vthr) < nanmean(spont.spkwave(:,1:end-1)),1,'first');
        axis([nanmean(spont.spk_Vthr)+[-2,2],tmpmean(a)+[-tmpmean(a),20]])
    end
    % Print figure
    
    set(fig,'paperposition',[0,0,18,15],'papersize',[18,15],'paperunits','centimeters')
    print(fig,'-dpdf',sprintf('sinusoidal_mod_%s.pdf',filename))
    close(fig)
    % Caption
    caption = sprintf(['Experiment %s. Sinusoidal modulation measured with a ',...
        ' %s waveform (nanmean %3.1f std %3.1f tau %3.1fms). Frequency of ',...
        'modulation %4.0f Hz, %3.1fpA amplitude and duration %3.1fs. ',...
        'The estimated capacitance from the noisy trace is: %4.1fpF. ',...
        'The spontaneous firing rate was %3.1fHz with a cv of %2.2f. ',...
        'Evoked firing rate was %3.1fHz and a total of %d spikes were',...
        ' used. Other parameters found on file.'],data.expName,...
        stim.type,stim.mean,stim.sigma,stim.tau,stim.F...
        ,stim.m,stim.tend-stim.tstart,spont.C,spont.frate,...
        spont.cv,r0,length(relspks));
    %
    printFigWithCaption(sprintf('sinusoidal_mod_%s.pdf',filename),caption)
    if exist(sprintf('sinusoidal_mod_%sCaption.pdf',filename),'file')
        movefile(sprintf('sinusoidal_mod_%sCaption.pdf',filename),sprintf('sinusoidal_mod_%s.pdf',filename))
    end
end % found_sinusoid

function [t, V, I, W, metadata, info, Vd, Id] = i_load_compensated_voltage(file,kfiles,holding,dendrite)
% Internal function to load voltage trace
V = [];
I = [];
Vd = [];
Id = [];
metadata = {};
[ent, info] = load_h5_trace(file.path);
idx = find(strcmp('RealNeuron',{ent.name}));
idx = [idx, find(strcmp('AnalogInput',{ent.name}))];


if ~isempty(dendrite)
    idx = idx(dendrite);
    V = [ent(idx(1)).data];
    Vd = [ent(idx(end)).data];
else
    
    V = [ent(idx(1)).data];
end

t = linspace(0,info.tend,length(V));

% Do we need AEC?
wave = find(strcmp('Waveform',{ent.name}));
W = ent(wave);
if isempty(ent(idx(1)).metadata)
    if length(dendrite) < 2
    % if there was a holding potential include it on the AEC current
        idx = [wave,find(strcmp('Constant',{ent.name}))];
        I = sum(vertcat(ent(idx).data),1) + holding(1);
        
        [~,k]  = min((file.date) - [kfiles.date]);
        Ke=load(kfiles(k).path);
        V = AECoffline(V,I,Ke);
        metadata = {ent(idx).metadata};
    else
        idx = find(strcmp('AnalogOutput',{ent.name})); % QUICKFIX!
        idx = idx(dendrite);
        I = ent(idx(1)).data + holding(1);
        Id = ent(idx(end)).data + holding(2); 
        metadata = {ent(idx).metadata};
        W = ent(idx(end));
        
        Ks = load(kfiles(dendrite(1)).path);
        Kd = load(kfiles(dendrite(2)).path);
        

        V = AECoffline(V,I,Ks);
        Vd = AECoffline(Vd,Id,Kd);
        
       
    end
 

end



function [spk,correction_value] = int_spikes...
    (spk,spikes,dt,int_factor,spike_alignment_mode)
% Interpolates spikes and returns corrected spike time.
%

nspk=size(spikes,1); % number of spikes
t=(0:size(spikes,2)-1).*dt; % samples space
intt=linspace(0,t(end),size(spikes,2)*int_factor); % interpolation space

intspikes=nan(nspk,length(intt));
correction_value=nan(size(spk));

for ii=1:nspk
    if sum(isnan(spikes(ii,:)))
        keyboard
    end
    intspikes(ii,:) = spline(t,spikes(ii,:),intt);
end
[~,maxspikes] = max(spikes,[],2);
[~,maxintspikes] = max(intspikes,[],2);
RMSFACTOR = 0.5;
switch spike_alignment_mode
    % interpolate spike shapes and correct/align spike times
    case 'thresh'
        %         d2Xi = diff(intspikes,2,2)./(dt*2);
        %       Uses the peak of the third derivative
        d3Xi = diff(intspikes,3,2)./(dt*3);
        for ii = 1:nspk
            [~,i3dV]        = findpeaks(d3Xi(ii,:),'minpeakdistance',3,'minpeakheight',rms(d3Xi(ii,:))*RMSFACTOR);
            %find the closes value to the most spike peak....but don't
            %correct for the number of samples (Trick!! As there is no true reason to do this??)
            if isempty(find(i3dV<maxintspikes(ii),1,'last'))
                fprintf(1,'There was an error in the spiketime correction!!!\n\n');
                keyboard
            else
                i3dV = i3dV(find(i3dV<maxintspikes(ii),1,'last'))-3;
            end
            %             plot(d3Xi(ii,:),'k')
            %             plot(i3dV,d3Xi(ii,i3dV),'ko')
            %             plot(xlim,rms(d3Xi(ii,:))*RMSFACTOR*[1,1],'r')
            %
            %             [~,i2dV]        = findpeaks(d2Xi(ii,:),'minpeakdistance',size(d2Xi,2)-2);%'minpeakheight',rms(d2Xi(ii,:))*RMSFACTOR,
            %             maxintspikes(ii) = i2dV(1)-2;
        end
end
% try
correction_value(:) = intt(maxintspikes)-t(maxspikes);
% catch
%     fprintf('There was an error in the threshold stuff...')
%     keyboard
% end
spk = spk + correction_value;