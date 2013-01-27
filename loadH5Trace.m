function [entities,info] = loadH5Trace(filename)
% [entities,info] = loadH5Trace(filename)
warning off
try
    version = hdf5read(filename, '/Info/version');
catch
    version = 1;
    info = hdf5info(filename);
    ngroups = length(info.GroupHierarchy.Groups);
    for k=1:ngroups
        if strcmp(info.GroupHierarchy.Groups(k).Name, '/Metadata');
            version = 0;
            break;
        end
    end
end

%fprintf(1, 'H5 file version #%.0f.\n', version);

switch version
    case 0
        [entities,info] = loadH5TraceV0(filename);
    case 1
        [entities,info] = loadH5TraceV1(filename);
    case 2
        [entities,info] = loadH5TraceV2(filename);
    otherwise
        error('Unknown H5 file version (%.0f).', version);
end
warning on
