function compute_prc_from_eLif(param,C,dt,Inoise,I0,Pamp,Pwidth)
Pamp = 150; %pA
Pwidth = 0.001; %s
F = [];
R = [];
figure(1),clf
ax(1) = axes('position',[0.1,0.1,0.5,.8]);
hold on
for ii = 1:length(I0)
    [x,y,frate] = i_compute_one_prc(param, C, dt, Inoise,I0(ii), Pamp, Pwidth);
    if frate~=0
        [ptb,me,ml,ie,il] = compute_peak_to_baseline(x, y);
        R = [R,ptb];
        F = [F,frate];
        
        plot(x,y.*1e3+.5*(ii-1),'k')
        plot(x(ie),y(ie).*1e3+.5*(ii-1),'ko')
        plot(x(il),y(il).*1e3+.5*(ii-1),'ko')
        plot([0,1],[0,0]+.5*(ii-1),'--','color',[0.5,.5,.5])
        text(0.1,.5*(ii),sprintf('%3.f Hz',frate),'color','r')
    end
end
set(gca,'ycolor','w','ytick',[])
xlabel('\Phi')
axis tight
xlim([0,1])
ax(2) = axes('position',[0.7,0.1,0.25,.35]);
plot(F,R,'ko')
xlabel('Frequency (Hz)')
ylabel('Peak to baseline')
ylim([0,1])
ax(3) = axes('position',[0.7,0.5,0.25,.35]);
[~,~,~,t,Vref,V] = i_compute_one_prc(param, C, dt, Inoise,I0(ceil(length(I0)/2)), Pamp, Pwidth);
plot(t*1e3,V,'r')
hold on
plot(t*1e3,Vref,'k')
axis tight
xlabel('time (ms)')
ylabel('V (mV')
set(gcf,'paperposition',[0,0,9,10],'papersize',[9,10],'paperunits','centimeters')
print(gcf,'-dpdf','eLIF_model_prc.pdf')


function [phi,dphi,frate,t,Vref,V] = i_compute_one_prc(param, C, dt, Inoise,I0, Pamp, Pwidth)
% computes the Phase Response Curve from an eLIF neuron.

Treset = 0.5*1e-3;
NPERT = 200;
% Run once. Does the model fire?
tmax = 1; %350 ms
phi = [];
dphi = [];
frate = 0;
t = 0:dt:tmax;

Vref = integrate_eLIF(param,t,I0*ones(size(t))+0*rand(size(t)),C,param(2),Treset);

spks = find(Vref>0);

if length(spks)<4
    disp('Error...Not enough spikes..')
    return
end

if length(spks)>8
    tmax = t(spks(8));
end
t = 0:dt:tmax;
Vref = Vref(find(t<=t(end)));
spks = find(Vref>0);
mean_isi = mean(diff(t(spks)*1e3));

p_size = Pwidth/dt;
min_stim_loc = spks(4);
max_stim_loc = spks(4)+mean(diff(spks));
pert_loc = unique(ceil(linspace(min_stim_loc,max_stim_loc,NPERT)));

Ptime = nan(size(pert_loc));
trials = cell(size(pert_loc));
spk_advance = nan(size(pert_loc));
tau = nan(size(pert_loc));

for ii = 1:length(pert_loc)
    I = I0*ones(size(t));
    I(pert_loc(ii):pert_loc(ii)+p_size) = Pamp+I0;
    V = integrate_eLIF(param,t,I+Inoise*(rand(size(t))-0.5),C,param(2),Treset);
    tmpspks = t(find(V>0))*1e3;
    Ptime(ii) = t(pert_loc(ii))*1e3;
    trials{ii} = tmpspks;
    spk_before = tmpspks(find(tmpspks<=Ptime(ii),1,'last'));
    spk_after = tmpspks(find(tmpspks>Ptime(ii),1,'first'));
    spk_advance(ii) = 1 - (spk_after - spk_before)/mean_isi;
    tau(ii) = Ptime(ii) - spk_before;
end
frate = 1000.0/mean_isi;
[tau,idx] = sort(tau);
phi = tau./mean_isi;
dphi = (spk_advance(idx))./(Pamp*Pwidth*1e3);

I = I0*ones(size(t));
ii = (ceil(NPERT/2));
I(pert_loc(ii):pert_loc(ii)+p_size) = Pamp+I0;
V = integrate_eLIF(param,t,I+Inoise*(rand(size(t))-0.5),C,param(2),Treset);

% ax1 = axes('position',[0.1,.5,.8,.4]);
% plot(t,Vref,'k-')
% hold on
% plot(t(spks),Vref(spks),'ko')


% plot(t,V,'r')
% axis tight
% plot(pert_loc(ii)*[1,1],ylim,'k--')
% ax2 = axes('position',[0.1,.1,.35,.35]);
% [tau,idx] = sort(tau);
% plot(tau,spk_advance(idx),'-')


