function [clu, tree,cluster_input] = run_superparamagnetic_clustering(inspk,cluster_input)
% [clu, tree,cluster_input] = run_superparamagnetic_clustering(inspk,cluster_input)

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
    cluster_input.tempstep = 0.001;                %temperature step
    cluster_input.num_temp = floor(...
        (cluster_input.maxtemp - ...
    cluster_input.mintemp)/cluster_input.tempstep); %total number of temperatures
    cluster_input.stab = 0.8;                     %stability condition for selecting the temperature
    cluster_input.SWCycles = 100;                 %number of montecarlo iterations
    cluster_input.KNearNeighb = 11;               %number of nearest neighbors
    cluster_input.randomseed = 0;                 % if 0, random seed is taken as the clock value
    %cluster_input.randomseed = 19;              % If not 0, random seed
    cluster_input.inputs = 10;
    dim=cluster_input.inputs;
end

tree = nan(cluster_input.num_temp,16);
clus = nan(cluster_input.num_temp,size(inspk,1)+2);
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
if exist([fname '.dg_01.lab'],'file')
    clu=load([fname '.dg_01.lab']);
    tree=load([fname '.dg_01']);
    delete(sprintf('%s.run',fname));
    delete *.mag
    delete *.edges
    delete *.param
    delete(fname_in);
end
