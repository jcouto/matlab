function varargout = loadH5Traces(files, idx)
% LOADH5TRACES
% [t , 1 data Array Per Entity , protocols] = loadH5Traces(files,entityIDs = [1,2], protocol)
% Returns a data array with the traces of the entities specified in idx.
% Each row is a data vector.
%
% idx are it indexes of the entities recorded. If a cell array of strings
% is used the strings should name the entities. Default is
% {'RealNeuron','Waveform'}.
% The last output contains the metadata from the specified entities.
%


% Default values
if ~exist('idx','var'),idx = {'RealNeuron','Waveform'};end

[entities, info] = loadH5Trace(files{1});

nent = length(idx);

% Locate entities names
if iscellstr(idx)
    tmp = idx;
    idx = zeros(length(idx),1);
    for ii = 1:length(tmp)
        index = find(~cellfun(@isempty,strfind({entities.name},tmp{ii})));
        if ~isempty(index)
            idx(ii) = index; 
        else
           warning('Entity %s not found (returning NaN).\n',tmp{ii})
           idx(ii)=0;
        end
    end
end

m = length(files);
n = length(entities(idx(1)).data);

% Generate time vector
t = (0:n-1) * info.dt;
metadata = cell(m,nent);

varargout = cell(1,nent+1);
% Create arrays to store the data
for ii = 1:nent
    varargout{1+ii} = nan(m,n);
    if idx(ii) > 0
        metadata{1,ii} = entities(idx(ii)).metadata;
        varargout{1+ii}(1,:) = entities(idx(ii)).data(:);
    end
end

for k=2:m
    [entities, ~] = loadH5Trace(files{k});
    for ii = 1:nent
        if idx(ii) > 0
            varargout{1+ii}(k,:) = entities(idx(ii)).data(:);
            metadata{k,ii} = entities(idx(ii)).metadata;
        end
    end
end
varargout{1}=t;
varargout{end+1}=metadata;


%Uncomment to display warning if the outputs don't match the requested in idx.
%if nargout<length(varargout)
%    warning('You specified too less output arrays. Use n_entities+2 as a thumb rule.')
%end

