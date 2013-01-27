function [entities,info] = loadH5TraceV1(filename)
% [entities,info] = loadH5TraceV1(filename)

info = hdf5info(filename);
ngroups = length(info.GroupHierarchy.Groups);

for ii=1:ngroups
    if strcmp(info.GroupHierarchy.Groups(ii).Name,'/Data') == 1
        nentities = length(info.GroupHierarchy.Groups(ii).Datasets);
        entities = repmat( ...
            struct('id',[],'data',[],'metadata',[],'units','','name',''), ...
            [nentities,1]);
        for jj=1:nentities
            name = info.GroupHierarchy.Groups(ii).Datasets(jj).Name;
            entities(jj).id = str2double(name(end-3:end));
            entities(jj).data = hdf5read(filename, name);
            nattrs = length(info.GroupHierarchy.Groups(ii).Datasets(jj).Attributes);
            for kk=1:nattrs
                name = info.GroupHierarchy.Groups(ii).Datasets(jj).Attributes(kk).Name;
                idx = strfind(name, '/');
                name = lower(name(idx(end)+1:end));
                if length(name) > 8 && strcmp(name(1:8),'metadata')
                    entities(jj).metadata = ...
                        info.GroupHierarchy.Groups(ii).Datasets(jj).Attributes(kk).Value';
                else
                    value = info.GroupHierarchy.Groups(ii).Datasets(jj).Attributes(kk).Value;
                    if isa(value, 'hdf5.h5string')
                        entities(jj).(name) = value.Data;
                    else
                        entities(jj).(name) = value;
                    end
                end
            end
        end
    end
end

clear('info');
info.tend = hdf5read(filename,'/Misc/Simulation_properties/tend');
info.dt = hdf5read(filename,'/Misc/Simulation_properties/dt');
info.srate = 1.0 / info.dt;
info.version = 1;
