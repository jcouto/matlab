function [class,tree,cluster_input,inspk, clu] = compute_spike_sorting_classes(spikes,...
    spike_wavelets,cluster_input, clus,tree,plotvar)
% [class,tree,cluster_input,inspk,clu] = compute_spike_sorting_classes(spikes,cluster_input,cluster,cluster_tree)
% The last class is the discarded waveforms.
%
%
% Adapted from WaveClus
if isempty(spikes)
    fprintf(1,'No spikes to cluster.\n');
    class = [];
    tree = [];
    cluster_input = [];
    inspk = [];
    return
end


if ~exist('spike_wavelets','var')
    inspk = [];
else
    inspk = spike_wavelets;
end
if isempty(inspk)
    [inspk] = compute_features_from_spk_shape(spikes);
end
if ~exist('clu','var') || ~exist('tree','var')
    clu = [];
    tree = [];
end
if isempty(clu) || isempty(tree)
    [clu,tree,cluster_input] = run_superparamagnetic_clustering(inspk);
end

[temp] = find_tree_temperature(tree,cluster_input,size(spikes,1),[],[],[]);

%%
% keyboard
%%
if exist('plotvar','var')
    n = size(clu,1);
    for t = 1:n
        labels = sort(unique(clu(t,3:end)));
        label_num = arrayfun(@(x)length(clu(t,(clu(t,3:end)==x))),labels);
        if max(label_num)<10
            break
        end
    end
    n = ceil(t);
    cc = setFigureDefaults();
    set(0,'DefaultFigureVisible','off'); 
    fig = figure('visible','off','papersize',[20,20],...
        'paperposition',[-1,-1,21,21],'paperunits','centimeter');
    boxon = 0;
    nplots_x = 4;
    nplots_y = (n/nplots_x)+1;
    count_x = 0;
    count_y = 0;
    for t = 1:n
        ax(1) = axes('position',[count_x./(nplots_x),...
            0.01+count_y./(nplots_y), 1./(nplots_x*2),...
            1./(nplots_y)]);
        ax(2) = axes('position',[count_x./(nplots_x)+1/(2*nplots_x),...
            0.01+count_y./(nplots_y) 1./(nplots_x*2),...
            1./(nplots_y)]);
        class = get_classes_for_temperature(clu,t);
        set(ax,'visible','off','box','off','drawmode','fast')
        plot_spike_shapes(class,spikes,inspk,ax,cc)
        
        axes(ax(1)), axis tight
        axes(ax(2)), axis tight
        tbox = annotation('textbox',[count_x./(nplots_x),...
                0.01+count_y./(nplots_y) 1./(nplots_x),...
                1./(nplots_y+1)]);
        set(tbox,'edgecolor','none','string',t,...
            'fontsize',12,'color','k','fontweight','bold'...
                ,'horizontalalignment','left','verticalalignment','bottom')
        if boxon
            boxon = 0;
            set(tbox,'edgecolor','k')
        else
            boxon = 1;
        end
        if t == temp
            set(tbox,'edgecolor','r','linewidth',2)
            boxon = 0;
        end
        count_x = count_x + 1;
        
        if count_x == nplots_x
            count_x = 0;
            count_y = count_y + 1;
        end
        clear('class','ax')
    end
    if ischar(plotvar)
        
        %print(fig,'-dpdf',sprintf('%s.pdf',plotvar))
        tbox = annotation('textbox',[0.1,0.9,0.8,0.1])
        set(tbox,'edgecolor','none','string',plotvar,...
            'fontsize',12,'color','k','fontweight','bold'...
                ,'horizontalalignment','center',...
                'verticalalignment','bottom',...
                'interpreter','none')
        print(fig,'-dpng',sprintf('%s.png',plotvar),'-r250')
        close(fig);
    end
end

class = get_classes_for_temperature(clu,temp);
%%
function plot_spike_shapes(class,spikes,inspk,ax,cc)
axes(ax(1))
plot(spikes(class{end},:)','-','linewidth',0.5,'color',[.5,.5,.5])
axes(ax(2))
plot(inspk(class{end},1),inspk(class{end},2),'ko',...
    'markeredgecolor',[.5,.5,.5],'markerfacecolor',[.5,.5,.5])
hold on
for ii = 1:length(class)-1
    axes(ax(1))
    plot(spikes(class{ii},:)','color',cc(ii,:),'linewidth',0.5)
    axes(ax(2))
    plot(inspk(class{ii},1),inspk(class{ii},2),'ko',...
        'markeredgecolor',cc(ii,:),'markerfacecolor','none')
end


function class = get_classes_for_temperature(clu,temp)
labels = sort(unique(clu(temp,3:end)));
label_num = arrayfun(@(x)length(clu(temp,(clu(temp,3:end)==x))),labels);
min_spk = 10;
labels = labels(label_num>min_spk);
class = cell(length(labels)+1,1);
selected_clu = clu(temp,3:end);
parfor ii = 1:length(labels)
    class{ii}=find(selected_clu==labels(ii));
    % plot(spikes(class{ii},:)','color',cc(ii),'linewidth',0.5)
end
% size(class)
class{end}= setdiff(1:size(clu,2)-2,[class{:}]);
% plot(spikes(class{end},:)','color',cc(ii),'linewidth',0.5)
% set(gcf,'visible','on')
