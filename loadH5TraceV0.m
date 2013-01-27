function [entities,info] = loadH5TraceV0(filename)
% out = loadH5TraceV0(filename)

info = hdf5info(filename);
ngroups = length(info.GroupHierarchy.Groups);

for ii=1:ngroups
    if strcmp(info.GroupHierarchy.Groups(ii).Name,'/Data') == 1
        nentities = length(info.GroupHierarchy.Groups(ii).Datasets);
        entities = repmat(struct('id',[],'data',[],'metadata',[],'parameters',[]), [nentities,1]);
        ids = zeros(nentities,1);
        for jj=1:nentities
            ids(jj) = str2double(info.GroupHierarchy.Groups(ii).Datasets(jj).Name(end-3:end));
            entities(jj).id = ids(jj);
            entities(jj).data = hdf5read(filename, info.GroupHierarchy.Groups(ii).Datasets(jj).Name);
        end
    elseif strcmp(info.GroupHierarchy.Groups(ii).Name,'/Metadata') == 1
        for jj=1:length(info.GroupHierarchy.Groups(ii).Datasets)
            name = info.GroupHierarchy.Groups(ii).Datasets(jj).Name;
            idx = strfind(name, '-');
            ntt = find(ids == str2double(name(idx(1)+1:idx(2)-1)));
            entities(ntt).metadata = hdf5read(filename, name)';
        end
    elseif strcmp(info.GroupHierarchy.Groups(ii).Name,'/Parameters') == 1
        for jj=1:length(info.GroupHierarchy.Groups(ii).Datasets)
            name = info.GroupHierarchy.Groups(ii).Datasets(jj).Name;
            ntt = find(ids == str2double(name(end-3:end)));
            entities(ntt).parameters = hdf5read(filename, name);
        end
    end
end

clear('info');
info.tend = hdf5read(filename,'/Misc/Simulation_properties/tend');
info.dt = hdf5read(filename,'/Misc/Simulation_properties/dt');
info.srate = 1.0 / info.dt;
info.version = 0;
