function process_sinusoidal_modulation_folder(folder)

if ~exist('folder','var')
    folder = pwd;
end
previous_folder = pwd;
cd(folder)
%% Compute the sin modulation if needed
files = list_h5_files();
data_files = list_files('.','sin_mod_data*',1);
data_files =[ data_files;list_files('.','no_modulation_data*',1)];

for ff = 1:length(files)
    [~,filename] = fileparts(files(ff).basename);
    if (sum(~cellfun(@isempty,strfind(data_files,filename)))<1)
        fprintf(1,'Processing file [%d/%d]: %s\n',ff,length(files),files(ff).basename);
        [~,found_sinusoid] = process_sinusoidal_modulation_file(files(ff).path);
        if ~found_sinusoid
            fprintf(1,'File [%d/%d]: %s .With did not have modulation.\n',ff,length(files),files(ff).basename);
        end
    end
end
%% Gather all
data_files = list_files('.','sin_mod_data*',1);
data = [];
for ff = 1:length(data_files)
    data = [data, load(data_files{ff})];
end
no_mod_data_files =list_files('.','no_modulation_data*',1);
no_mod_data = [];
for ff = 1:length(no_mod_data_files)
    no_mod_data = [no_mod_data, load(no_mod_data_files{ff})];
end
%% Extract frequency and spontaneous parameters

% Modulation
r0 = [data.r0];
r1 = [data.r1];
phi = [data.phi];
conf = [data.conf];
shuffled_r0 = [data.r0_shuffled];
shuffled_r1 = [data.r1_shuffled];
shuffled_phi = [data.phi_shuffled];
r = [data.r];
N = nan(length(data),1);
F = nan(length(data),1);
% Capacitance estimate
C = nan(length(data),1);
% spontaneous spikes
spk.frate = nan(length(data),1);
spk.height = nan(length(data),1);
spk.width = nan(length(data),1);
spk.slope = nan(length(data),1);
spk.thresh = nan(length(data),1);
spk.ahp = nan(length(data),1);
% Stimulation ASSUMES THEY ARE ALL THE SAME
stim.sigma = data(1).stim.sigma;
stim.mean = data(1).stim.mean;
stim.tau = data(1).stim.tau;
% Coefficient of variation (in the absence of stimulation)
for ii = 1:length(data)
    F(ii) = data(ii).stim.F;
    N(ii) = length(data(ii).spiketimes(...
        data(ii).spiketimes>data(ii).stim.tstart &...
        data(ii).spiketimes<data(ii).stim.tend));
    C(ii) = data(ii).spont.C;
    spk.height(ii) = mean(data(ii).spont.spk_height);
    spk.width(ii) = mean(data(ii).spont.spk_width);
    spk.slope(ii) = mean(data(ii).spont.spk_slope);
    spk.thresh(ii) = mean(data(ii).spont.spk_Vthr);
    spk.ahdp(ii) = mean(data(ii).spont.spk_ahdp);
    spk.frate(ii) = nanmean(1./diff(data(ii).spont.spks));
end
% Sort indexes for the frequency
[~,idx] = sort(F);
% From no_mod_data
if ~isempty(no_mod_data)
coeffVar = no_mod_data.cv;
firingRate = no_mod_data.fr;
capacitance = no_mod_data.spont.C;
else
    coeffVar = nan;
    firingRate = nan;
    capacitance = nan;
end
expName = data(1).expName;
%% Plot figure with all the histograms on the same page.
cc = setFigureDefaults;
ax = [];
counter = 1;
nBins = size(r,1);
figCounter = 0;
for ii = 1:length(F)
    if counter == 1
        fig = figure();
        figCounter = figCounter+1;
    end
    ax = [ax, subplot(4,3,counter)];
    edges = linspace(0,1000./F(idx(ii)),length(r(:,(idx(ii)))));
    x = linspace(0,1000./F(idx(ii)),nBins*10);
    hold on;
    plot([edges(1),edges(end)],r0(idx(ii))+[0,0],'--','color',cc(2,:),'LineWidth',1.2);
    
    hndl = bar(edges, r(:,(idx(ii))), 0.8, 'histc');
    
    set(hndl,'FaceColor',[0.5,.5,.5],'linewidth',1);
    plot(x, r0(idx(ii))+r1(idx(ii))*sin(2*pi*F(idx(ii))*x/1000.0), '--','color',cc(3,:), 'LineWidth', 1.5);
    plot(x, r0(idx(ii))+r1(idx(ii))*sin(2*pi*F(idx(ii))*x/1000.0+phi(idx(ii))),'color',cc(1,:), 'LineWidth', 1.2);
    %     axis([0,edges(end),min(r(1:end-1,(idx(ii))))-5,max(r(:,(idx(ii))))+5]);
    axis('tight')
    set(ax,'ylim',[0,max(r(:))+5])
    text(min(ylim),max(ylim),sprintf('\t\tF = %4dHz; N = %d spikes.',F(idx(ii)),N(idx(ii))),'verticalalignment','top','horizontalalignment','left')
    counter = counter + 1;
    if counter == 13 | ii == length(F)
        xlabel({'Time from sinusoidal cycle (ms)'})
        ylabel('Firing frequency (Hz)')
        counter = 1;
        filename = sprintf('sinusoidal_modulation_%s_%02d',expName,figCounter);
        set(fig,'paperposition',[0,0,18,15],'papersize',[18,15],'paperunits','centimeters')
        print(fig,'-dpdf',sprintf('%s.pdf',filename))
        figCounter = figCounter+1;
        caption = sprintf('Experiment %s. Histogram of the sinusoidal modulation for frequencies ',data(1).expName);
        for jj = 1:length(F)-1
            caption = sprintf('%s %d,',caption,F(jj));
        end
        caption = sprintf('%s and %d Hz. The red line is the fit and the blue dashed the stimulus.',caption,F(end));
        printFigWithCaption(sprintf('%s.pdf',filename),caption)
        if exist(sprintf('%sCaption.pdf',filename),'file');
            movefile(sprintf('%sCaption.pdf',filename),sprintf('%s.pdf',filename))
        end
        close(fig)
    end
end

%% Compute the Transfer Function fit
% adjust the values of phi be monotonically decreasing
tmpphi = phi(idx);
start = find(tmpphi < 0, 1);
for k=start+1:length(tmpphi)
    while tmpphi(k) > tmpphi(k-1)+pi/2
        tmpphi(k) = tmpphi(k) - 2*pi;
    end
end
if tmpphi(1) < -pi
    tmpphi = tmpphi+2*pi;
end

w = [diff(reshape([conf.r1],[2,length([conf.r1])/2]))./diff(reshape([conf.r0],[2,length([conf.r0])/2]));diff(reshape([conf.phi],[2,length([conf.phi])/2]))];
[A,Z,P,Dt,Mopt,Nopt,fvar] = fitTransferFunction(F(idx), (r1(idx)./r0(idx))',tmpphi',[1:2], [2:3],w(:,idx)');


%% save data

savedata = spk;
savedata.expName = expName;
% Modulation
savedata.r0 = r0;
savedata.r1 = r1;
savedata.phi = phi;
savedata.corrected_phi = ipermute(tmpphi,idx);
savedata.conf = conf;
savedata.shuffled_r0 = shuffled_r0;
savedata.shuffled_r1 = shuffled_r1;
savedata.shuffled_phi = shuffled_phi;
savedata.r = r;
savedata.N = N;
savedata.F = F;
savedata.C = C;
savedata.sigma = stim.sigma;
savedata.mean = stim.mean;
savedata.tau = stim.tau;
savedata.sort_F_idx = idx;
savedata.CV = coeffVar; 
savedata.frate = firingRate;
savedata.capacitance = capacitance;
savedata.A = A;
savedata.P = P;
savedata.Z = P;
savedata.Dt = Dt;
savedata.fit_eval = fvar;

save(sprintf('processed_modulation_data%s.mat',expName),'-struct','savedata');

%%
fig = figure();
ax(1) = axes('position',[0.1,0.4,.4,.25],'yscale','linear','xscale','log');
plot(fvar.f,10*log10(fvar.mag),'-','color',cc(1,:),'LineWidth',1.2);
hold on;
plot(F(idx),10*log10(r1(idx)./r0(idx)),'ko','MarkerFaceColor','k');
plot(F(idx),10*log10(shuffled_r1(idx)./shuffled_r0(idx)),'o',...
    'Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5]);
axis tight;
for i=1:length(Z)
    plot(abs(Z(i))+[0,0],ylim,'--','color',cc(2,:));
end
for i=1:length(P)
    plot(abs(P(i))+[0,0],ylim,'--','color',cc(3,:));
end
if ~isempty(Z(find((Z>1 & Z<2000),1,'last')))
text(Z(find((Z>1 & Z<2000),1,'last')),min(ylim),sprintf('\t\tZeros (%d)',Mopt),'color',cc(2,:),...
    'verticalalignment','bottom','horizontalalignment',...
    'left','rotation',90);
end
if ~isempty(P(find((P>1 & P<2000),1,'last')))
text(P(find((P>1 & P<2000),1,'last')),min(ylim),sprintf('\t\tPoles (%d)',Nopt),'color',cc(3,:),...
    'verticalalignment','bottom','horizontalalignment',...
    'left','rotation',90);
end
xlabel('Frequency (Hz)')
ylabel('10log_{10}({r_1/r_0}) (dB)');

ax(2) = axes('position',[0.1,0.7,.4,.25],'yscale','linear','xscale','log',...
    'ycolor','k','yaxislocation','left');

plot(fvar.f,abs(fvar.mag),'-','color',cc(4,:),'LineWidth',1.2);
hold on;
%
plot(F(idx),(r1(idx)./r0(idx)),'x','MarkerEdgeColor','k');
plot(F(idx),(shuffled_r1(idx)./shuffled_r0(idx)),'x',...
    'Color',[.5,.5,.5],'MarkerEdgeColor',[.5,.5,.5]);
axis tight;
ylabel('(r_1/r_0)');
for i=1:length(Z)
    plot(abs(Z(i))+[0,0],ylim,'--','color',cc(2,:));
end
for i=1:length(P)
    plot(abs(P(i))+[0,0],ylim,'--','color',cc(3,:));
end
ylabel('(r_1/r_0)');

ax(3) = axes('position',[0.1,0.1,.4,.2],'yscale','linear','xscale','log');

% Plot phase
plot(F(idx),tmpphi*180/pi,'ko-','LineWidth',1,'MarkerFaceColor','k');
hold on;
plot(F(idx),180/pi*shuffled_phi(idx),'bo','MarkerFaceColor',[.5,.5,.5]);
plot(fvar.f,fvar.phi*180/pi,'-','color',cc(1,:),'LineWidth',1.2);
xlabel('Frequency (Hz)')
ylabel('Phase(\Phi)');
for i=1:length(Z)
    plot(abs(Z(i))+[0,0],ylim,'--','color',cc(2,:));
end
for i=1:length(P)
    plot(abs(P(i))+[0,0],ylim,'--','color',cc(3,:));
end
warning off
linkaxes(ax,'x')
warning on

ax(4) = axes('position',[0.6,0.4,.3,.2],'xtick',[]);
plot(C,'ko--','color','k','markerfacecolor','k');
ylabel('Capacitance estimate (pF)')
ylim([0,800])
ax(5) = axes('position',[0.6,0.4,.3,.2],'yaxislocation','right','ycolor',cc(1,:));
ylabel('Voltage (mV)')
plot(spk.thresh(idx),'ko--','color',cc(1,:),'markerfacecolor',cc(1,:));
text(max(xlim),max(spk.thresh),'threshold     ','color',cc(1,:),...
    'verticalalignment','bottom','horizontalalignment','right')
plot(spk.height(idx),'ko--','color',cc(2,:),'markerfacecolor',cc(2,:));
text(min(xlim),max(spk.height),'    height','color',cc(2,:),...
    'verticalalignment','bottom','horizontalalignment','left')


ax(6) = axes('position',[0.6,0.7,.3,.2]);
plot(spk.slope(idx),'ko--','color','k','markerfacecolor','k');
ylabel('Velocity (mV.ms-1)')

ax(7) = axes('position',[0.6,0.7,.3,.2],'yaxislocation','right','ycolor',cc(1,:));
plot(spk.width(idx),'ko--','color',cc(1,:),'markerfacecolor',cc(1,:));
ylabel('Spike width (ms)')
xlabel('Trial sorted by Frequency')

ylim([0.1,2])


ax(8) = axes('position',[0.6,0.1,.3,.2],'yaxislocation','left','ycolor','k');
plot(spk.frate(idx),'ko--','color','k','markerfacecolor',cc(1,:));
ylabel('Spontaneous Firing rate (Hz)')
xlabel('Trial sorted by Frequency')
ylim([0.1,1.5]*max(spk.frate))
linkaxes(ax([4:8]),'x')
xlim([0,length(idx)])

filename = sprintf('sinusoidal_modulation_%s_%02d',expName,figCounter);
set(fig,'paperposition',[0,0,18,15],'papersize',[18,15],'paperunits','centimeters')
print(fig,'-dpdf',sprintf('%s.pdf',filename))

caption = sprintf(['Experiment %s. Transfer function fit of the input output',...
    ' relation. Noisy stimulation with std %3.2f pA and tau %3.2f ms. Frequency',...
    ' modulated from %.0f to %4.0fHz. Fitted with rational transfer function (%d x %d), the constant',...
    ' term is %3.1f, the zeros'],...
    expName,stim.sigma,stim.tau,min(F),max(F),length(Z),length(P),A);
for zz = 1:length(Z)
caption = sprintf('%s %.2f,',caption,Z(zz));
end
caption = sprintf('%s and the poles ',caption);
for zz = 1:length(P)-1
caption = sprintf('%s %.2f,',caption,P(zz));
end
caption = sprintf('%s %.2f.',caption,P(end));
 caption = sprintf(['%s The CV of the data with no ',...
    'sinusoid is %.2f at %.1f Hz. The estimated capa',...
    'citance without modulation is %.1f pF. The ',...
    'figures to the right illustrate the evolution of',...
    ' the spontaneous spike parameters accross the ',...
    'experiment.'],caption,coeffVar,firingRate,capacitance);

printFigWithCaption(sprintf('%s.pdf',filename),caption)
if exist(sprintf('%sCaption.pdf',filename),'file');
    movefile(sprintf('%sCaption.pdf',filename),sprintf('%s.pdf',filename))
end

close(fig)

%%
cd(previous_folder)