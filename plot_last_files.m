function plot_last_files(N,mode)
% Plots the last N h5 files. 
% N is the number of files to plot (default N = 1) 
%
% Use the second input to specify the type of plot.
%   - voltage
%   - stack
% Joao Couto, May 2013

if ~exist('mode','var'); mode = 'voltage';end

[files,kfiles]=list_h5_files;
if(length(files) < N);N = length(files);end
switch mode
    case 'voltage'
        clf;
        ax(1) = axes();% Vm
        ax(2) = axes();% dVdt
        ax(3) = axes();% Vm hist
        ax(4) = axes();% IFF hist
        for ii = 1:N 
           [t,V,info] = i_load_compensated_voltage(files(end-1),kfiles);
           dVdt = diff(V)./(info.dt*1.e3);
           Vedges = min(V)-5:0.5:max(V)+5;
           Vbins = histc(V,Vedges);
           spk_idx = argfindpeaks(V,-20);
           
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
        set(ax,'box','off','tickdir','out')
        set(ax(1),'position',[0.1,0.1,0.55,0.4])
        set(ax(2),'position',[0.1,0.55,0.55,0.4],'xaxislocation','top',...
            'xcolor','w')
        set(ax(3),'position',[0.7,0.1,0.2,0.4])
        set(ax(4),'position',[0.75,0.6,0.2,0.3])
        set(ax([3,4]),'yaxislocation','right','xdir','reverse')
        linkaxes(ax([1,3]),'y');linkaxes(ax([1,2]),'x');
    case 'stack'
        'hello'
    otherwise
        print('Unknown mode')
end


function [t, V,info] = i_load_compensated_voltage(file,kfiles)
% Internal function to load voltage trace
    [ent, info] = loadH5Trace(file.path);
    idx = find(strcmp('RealNeuron',{ent.name}));
    idx = [idx, find(strcmp('AnalogInput',{ent.name}))]
    V = [ent(idx(1)).data];
    % Do we need AEC?
    if isempty(ent(idx(1)).metadata) 
        [~,k]  = min((file.date) - [kfiles.date]);     
        Ke=load(kfiles(k).path);
        V = AECoffline(V,I,Ke);
    end
    t = linspace(0,info.tend,length(V));