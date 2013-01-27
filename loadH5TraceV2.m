function [entities,info] = loadH5TraceV2(filename)
% [entities,info] = loadH5TraceV2(filename)

try
    version = hdf5read(filename, '/Info/version');
    if version ~= 2
        error('Unknown version in H5 file.');
    end
catch
    error('Unknown version in H5 file.');
end

info = hdf5info(filename);
ngroups = length(info.GroupHierarchy.Groups);

for ii=1:ngroups
    if strcmp(info.GroupHierarchy.Groups(ii).Name,'/Entities') == 1
        nentities = length(info.GroupHierarchy.Groups(ii).Groups);
        entities = repmat( ...
            struct('id',[],'data',[],'metadata',[],'units','','name','','parameters',struct([])), ...
            [nentities,1]);
        for jj=1:nentities
            name = info.GroupHierarchy.Groups(ii).Groups(jj).Name;
            idx = strfind(name, '/');
            entities(jj).id = str2double(name(idx(end)+1:end));
            entities(jj).data = hdf5read(filename, [name,'/Data']);
            try
                entities(jj).metadata = hdf5read(filename, [name,'/Metadata'])';
            catch
            end
            nattrs = length(info.GroupHierarchy.Groups(ii).Groups(jj).Attributes);
            for kk=1:nattrs
                name = info.GroupHierarchy.Groups(ii).Groups(jj).Attributes(kk).Name;
                idx = strfind(name, '/');
                name = lower(name(idx(end)+1:end));
                value = info.GroupHierarchy.Groups(ii).Groups(jj).Attributes(kk).Value;
                if isa(value, 'hdf5.h5string')
                    value = value.Data;
                end
                entities(jj).(name) = value;
            end
            nsubgroups = length(info.GroupHierarchy.Groups(ii).Groups(jj).Groups);
            for kk=1:nsubgroups
                name = info.GroupHierarchy.Groups(ii).Groups(jj).Groups(kk).Name;
                idx = strfind(name, '/');
                name = name(idx(end)+1:end);
                if strcmpi(name, 'parameters')
                    npars = length(info.GroupHierarchy.Groups(ii).Groups(jj).Groups(kk).Attributes);
                    for ll=1:npars
                        name = info.GroupHierarchy.Groups(ii).Groups(jj).Groups(kk).Attributes(ll).Name;
                        idx = strfind(name, '/');
                        name = lower(name(idx(end)+1:end));
                        value = info.GroupHierarchy.Groups(ii).Groups(jj).Groups(kk).Attributes(ll).Value;
                        if isa(value, 'hdf5.h5string')
                            value = value.Data;
                        end
                        entities(jj).parameters(1).(name) = value;                            
                    end
                end
            end
        end
    end
end

clear('info');
info.tend = hdf5read(filename,'/Info/tend');
info.dt = hdf5read(filename,'/Info/dt');
info.srate = 1.0 / info.dt;
info.version = 2;
