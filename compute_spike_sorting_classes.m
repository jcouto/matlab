function [class,cluster_input] = compute_spike_sorting_classes(spikes,cluster_input)
% Adapted from WaveClus

min_clus_abs = 60;              %minimum cluster size (absolute value)
min_clus_rel = 0.005;           %minimum cluster size (relative to the total nr. of spikes)
min_spks = 10;
[inspk] = compute_wavelets_from_spk_shape(spikes);
[clu,tree,cluster_input] = run_cluster(inspk);

[temp] = find_temp(tree,1,max(min_clus_abs,min_clus_rel*size(spikes,1)));
class = cell(5,1);
% cc = 'krbygb';
% hold all
for ii = 1:length(class)
    class{ii}=find(clu(temp,3:end)==ii-1);
%     plot(spikes(class{ii},:)','color',cc(ii))
end
% size(class)
class{ii+1}=setdiff(1:size(spikes,1), sort([class{:}]));

function [clu, tree,cluster_input] = run_cluster(inspk,cluster_input)

fname='tmp_SPC';
fname_in='spike_wavelets';
save(fname_in,'inspk','-ascii');
% DELETE PREVIOUS FILES
warning off
eval( ['delete([fname ''.dg_01.lab''])' ],'[];');
eval( ['delete([fname ''.dg_01''])' ],'[];');
warning on
if ~exist('cluster_input','var')
    cluster_input.mintemp = 0;                    %minimum temperature
    cluster_input.maxtemp = 0.201;                %maximum temperature
    cluster_input.tempstep = 0.01;                %temperature step
    cluster_input.num_temp = floor(...
        (cluster_input.maxtemp - ...
    cluster_input.mintemp)/cluster_input.tempstep); %total number of temperatures
    cluster_input.stab = 0.8;                     %stability condition for selecting the temperature
    cluster_input.SWCycles = 100;                 %number of montecarlo iterations
    cluster_input.KNearNeighb = 11;               %number of nearest neighbors
    %cluster_input.randomseed = 0;                 % if 0, random seed is taken as the clock value
    cluster_input.randomseed = 19;              % If not 0, random seed
    cluster_input.inputs = 10;
    dim=cluster_input.inputs;
end

dat=load(fname_in);
n=length(dat);
fid=fopen(sprintf('%s.run',fname),'wt');
fprintf(fid,'NumberOfPoints: %s\n',num2str(n));
fprintf(fid,'DataFile: %s\n',fname_in);
fprintf(fid,'OutFile: %s\n',fname);
fprintf(fid,'Dimensions: %s\n',num2str(dim));
fprintf(fid,'MinTemp: %s\n',num2str(cluster_input.mintemp));
fprintf(fid,'MaxTemp: %s\n',num2str(cluster_input.maxtemp));
fprintf(fid,'TempStep: %s\n',num2str(cluster_input.tempstep));
fprintf(fid,'SWCycles: %s\n',num2str(cluster_input.SWCycles));
fprintf(fid,'KNearestNeighbours: %s\n',num2str(cluster_input.KNearNeighb));
fprintf(fid,'MSTree|\n');
fprintf(fid,'DirectedGrowth|\n');
fprintf(fid,'SaveSuscept|\n');
fprintf(fid,'WriteLables|\n');
fprintf(fid,'WriteCorFile~\n');
if cluster_input.randomseed ~= 0
    fprintf(fid,'ForceRandomSeed: %s\n',num2str(cluster_input.randomseed));
end
fclose(fid);

[str,maxsize,endian]=computer;
cluster_input.system=str;
switch cluster_input.system
    case 'PCWIN'
        progpath = which('cluster_win.exe');
        runcmd = sprintf('%s %s.run',progname,fname);
    case 'MAC'
        progpath = which('cluster_mac.bin');
        runcmd = sprintf('%s %s.run',progpath,fname);
     case 'MACI64'
        progpath = which('cluster_maci.bin');
        runcmd = sprintf('%s %s.run',progpath,fname);
    otherwise  %(GLNX86, GLNXA64, GLNXI64 correspond to linux)
        progpath = which('cluster_linux.bin');
        run_cmd = sprintf('%s %s.run',progpath,fname);
end
[status,result] = system(runcmd);
if status
    disp('SPC failed...')
end
clu=load([fname '.dg_01.lab']);
tree=load([fname '.dg_01']);
delete(sprintf('%s.run',fname));
delete *.mag
delete *.edges
delete *.param
delete(fname_in);

function [temp] = find_temp(tree,num_temp,min_clus);
% Selects the temperature.

num_temp=num_temp;
min_clus=min_clus;

aux =diff(tree(:,5));   % Changes in the first cluster size
aux1=diff(tree(:,6));   % Changes in the second cluster size
aux2=diff(tree(:,7));   % Changes in the third cluster size
aux3=diff(tree(:,8));   % Changes in the third cluster size

temp = 1;         % Initial value

for t=1:num_temp-1;
    % Looks for changes in the cluster size of any cluster larger than min_clus.
    if ( aux(t) > min_clus | aux1(t) > min_clus | aux2(t) > min_clus | aux3(t) >min_clus )    
        temp=t+1;         
    end
end

%In case the second cluster is too small, then raise the temperature a little bit 
if (temp == 1 & tree(temp,6) < min_clus)
    temp = 2;
end   