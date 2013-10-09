function [temp_idx] = find_tree_temperature(tree, cluster_input, nspikes,min_clus_abs,min_clus_rel, min_clus,plotvar)
%[temp_idx] = find_tree_temperature(tree, cluster_input, nspikes,min_clus_abs,min_clus_rel, min_clus,plotvar)
% Selects the temperature of a spike sorting tree.
% The arguments are:
%    - tree: the cluster tree as computed by the superparamagnetic
%    clustering
%    - cluster_input: the input from the SPC parameters
%    - nspikes: number of spikes
%    - min_clus_rel: minimum cluster size (relative to the total nr. of spikes)
%    - min_clus_abs: minimum cluster size (absolute)
%    - specify the value of min_clus (overides the min_clus_rel and min_clus_rel)
%    - plotvar: Plots the temperature.
% Returns the temperature

if ~exist('min_clus_abs','var') || isempty('min_clus_abs') 
    min_clus_abs = 80;
end
if isempty(min_clus_abs)
    min_clus_abs = 80;
end
if ~exist('nspikes','var') || isempty('nspikes')
    nspikes = 6000;
end

if ~exist('min_clus_rel','var')
    min_clus_rel = 0.005;
end
if isempty(min_clus_rel)
    min_clus_rel = 0.005;
end

if ~exist('min_clus','var') 
    min_clus = max([min_clus_abs,min_clus_rel*nspikes]);
end

if isempty(min_clus)
    min_clus = max([min_clus_abs,min_clus_rel*nspikes])
end

num_temp = cluster_input.num_temp;

aux =diff(tree(:,5));   % Changes in the first cluster size
aux1=diff(tree(:,6));   % Changes in the second cluster size
aux2=diff(tree(:,7));   % Changes in the third cluster size
aux3=diff(tree(:,8));   % Changes in the third cluster size

temp_idx = 1;         % Initial value

for t=1:num_temp-1;
    % Looks for changes in the cluster size of any cluster larger than min_clus.
    if ( aux(t) > min_clus | aux1(t) > min_clus | aux2(t) > min_clus | aux3(t) > min_clus )    
        temp_idx=t+1;         
    end
end

%In case the second cluster is too small, then raise the temperature a little bit 
if (temp_idx == 1 & tree(temp_idx,6) < min_clus)
    temp_idx = 2;
end

if exist('plotvar','var')
    temperature=cluster_input.mintemp+temp_idx*cluster_input.tempstep;
    hold all
    plot([cluster_input.mintemp cluster_input.maxtemp-cluster_input.tempstep],[min_clus min_clus],'k:')
    plot(cluster_input.mintemp+(1:cluster_input.num_temp)*cluster_input.tempstep, ...
            tree(1:cluster_input.num_temp,5:size(tree,2)))
    plot([temperature temperature],[1 tree(1,5)],'k:')
    
end
