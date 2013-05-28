function [ent,info] = open_h5_trace(filename)
% This function opens a h5 file and returns pointers to the data.
% 
%

try
    key = '/Info';
    tmp = h5info(filename, key);
catch
    disp('Unknown format...')
    ent = [];
    info = [];
    return
end

% Extract "info" from attributes
for ii = 1:length(tmp.Attributes)
    info.(tmp.Attributes(ii).Name) = tmp.Attributes(ii).Value;
end
% And now the data and metadata
key = '/Entities';
tmp = h5info(filename, key);
for ii = 1:length(tmp.Groups)
    
    for jj = 1:length(tmp.Groups(ii).Attributes)
        ent(ii).(lower(tmp.Groups(ii).Attributes(jj).Name)) = tmp.Groups(ii).Attributes(jj).Value;
    end
    for jj = 1:length(tmp.Groups(ii).Datasets)
        ds_path = sprintf('%s/%s',tmp.Groups(ii).Name,tmp.Groups(ii).Datasets(jj).Name);
        ent(ii).(lower(tmp.Groups(ii).Datasets(jj).Name)) = h5trace_dataclass(filename,ds_path);
    end
end
end
