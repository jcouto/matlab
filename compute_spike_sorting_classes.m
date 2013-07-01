function [class,tree,cluster_input,inspk] = compute_spike_sorting_classes(spikes,spike_wavelets,cluster_input, cluster,cluster_tree)
% [class,tree,cluster_input,inspk] = compute_spike_sorting_classes(spikes,cluster_input,cluster,cluster_tree)
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
    [inspk] = compute_wavelets_from_spk_shape(spikes);
end
if ~exist('cluster','var') || ~exist('cluster_tree','var')
    clu = [];
    tree = [];
end
if isempty(clu) || isempty(tree)
    [clu,tree,cluster_input] = run_superparamagnetic_clustering(inspk);
    
end

[temp] = find_tree_temperature(tree,cluster_input,[],[],[],[]);
class = cell(5,1);
% cc = 'krbygb';
% hold all
for ii = 1:length(class)
    class{ii}=find(clu(temp,3:end)==ii-1);
%     plot(spikes(class{ii},:)','color',cc(ii))
end
% size(class)
class{ii+1}=setdiff(1:size(spikes,1), sort([class{:}]));
