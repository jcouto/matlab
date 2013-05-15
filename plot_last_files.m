function plot_last_files(N,mode,optvar)
% Plots the last N h5 files.
% N is the number of files to plot (default N = 1)
%
% Use the second input to specify the type of plot.
%   - voltage
%   - kink
%   - stack
% Joao Couto, May 2013

if ~exist('mode','var'); mode = 'voltage';end
if ~isstr(mode); optvar = mode;mode = 'voltage';end

[files,kfiles]=list_h5_files;
if(length(files) < N);N = length(files),end
cc = [21,  69, 152;...
    230,  20,  36;...
    19, 124,  56;...
    3,   4,   4;...
    239, 104,  28;...
    82,  24, 126;...
    143,  14,  25;...
    163,  31, 129;...
    96,  96,  96;...
    234,  66,  78;...
    196, 186,  87;...
    74, 136, 252;...
    247, 150,  73;...
    140,  79, 155;...
    140,  79, 155;...
    193,  90,  71;...
    204, 103, 165;...
    165, 165, 165]./256;
Ncc = size(cc,1);
switch mode
    case 'voltage'
        clf;
        ax(1) = axes();% Vm
        ax(2) = axes();% dVdt
        ax(3) = axes();% Vm hist
        ax(4) = axes();% IFF hist
        for ii = N:-1:1
            [t,V,~,info] = i_load_compensated_voltage(files(end-ii+1),kfiles);
            if ~exist('optvar','var')
                optvar = [0,t(end)];
            elseif length(optvar)<2
                optvar = [0,t(end)];
            end
            idx = (t>=optvar(1) &t<=optvar(2));
            t = t(idx);V = V(idx);
            dVdt = diff(V)./(info.dt*1.e3);
            Vedges = min(V)-5:0.5:max(V)+5;
            Vbins = histc(V,Vedges);
            spk_idx = argfindpeaks(V,-20);
            figure(1)
            axes(ax(1))
            plot(t,V,'k','linewidth',0.7)
            hold on
            plot(t(spk_idx),V(spk_idx),'ok','markerfacecolor','r')
            axes(ax(2))
            plot(t(2:end),dVdt,'k','linewidth',0.7)
            hold on
            plot(t(spk_idx-1),dVdt(spk_idx-1),'ok','markerfacecolor','r')
            axes(ax(3))
            barh(Vedges,Vbins,'k')
        end
        set(ax(1),'position',[0.1,0.1,0.55,0.4])
        set(ax(2),'position',[0.1,0.55,0.55,0.4],'xaxislocation','top',...
            'xcolor','w')
        set(ax(3),'position',[0.7,0.1,0.2,0.4])
        set(ax(4),'position',[0.75,0.6,0.2,0.3])
        set(ax([3,4]),'yaxislocation','right','xdir','reverse')
        linkaxes(ax([1,3]),'y');linkaxes(ax([1,2]),'x');
    case 'stack'
        'hello'
    case 'kink'
        clf;
        ax(1) = axes();% Vm
        ax(2) = axes();% dVdt
        for ii = N:-1:1
            [t,V,~,info] = i_load_compensated_voltage(files(end-ii+1),kfiles);
            if ~exist('optvar','var')
                optvar = [0,t(end)];
            elseif length(optvar)<2
                optvar = [0,t(end)];
            end
            tmin = optvar(1);tmax = optvar(2);
%             % find spk threshold
%             dVdt = diff(V)./(info.dt*1.e3);
%             idx = find(dVdt>100);
%             threshold = median(V(idx));
            [spk, spk_w, tspk_w] = extract_spikes( V(t>tmin&t<tmax), [], t(t>tmin&t<tmax), 4, 20, 20);
            dspk_w = diff(spk_w, 1, 2)./(info.dt*1e3);
            axes(ax(1))
            plot(spk_w(:,1:end-1)',dspk_w','color',cc(mod(ii,Ncc),:))
            hold on;
            axes(ax(2))
            plot(tspk_w,spk_w,'color',cc(mod(ii,Ncc),:))
            hold on;
        end
        set(ax(1),'position',[0.55,0.1,0.4,0.8])
        axes(ax(1));xlabel('time (ms)');ylabel('dV/dt (mV.ms^{-1})')
        set(ax(2),'position',[0.1,0.1,0.4,0.8])
        axes(ax(2));xlabel('time (ms)');ylabel('dV (mV)')
    otherwise
        print('Unknown mode')
end


function [t, V, I, info] = i_load_compensated_voltage(file,kfiles)
% Internal function to load voltage trace
V = [];
I = [];
[ent, info] = load_h5_trace(file.path);
idx = find(strcmp('RealNeuron',{ent.name}));
idx = [idx, find(strcmp('AnalogInput',{ent.name}))];
V = [ent(idx(1)).data];
t = linspace(0,info.tend,length(V));
% Do we need AEC?
if isempty(ent(idx(1)).metadata) && ~isempty(kfiles)
    % Is there I?
    idx = find(strcmp('Waveform',{ent.name}));
    idx = [idx,find(strcmp('Constant',{ent.name}))];
    I = sum(vertcat(ent(idx).data),1);
    [~,k]  = min((file.date) - [kfiles.date]);
    Ke=load(kfiles(k).path);
    V = AECoffline(V,I,Ke);
end
