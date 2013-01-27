function varargout = loadH5Traces(files,idx)
% LOADH5TRACES
% [t , 1 data Array Per Entity , protocols] = loadH5Traces(files,entityIDs = [1,2])

if isstruct(files(1))
    files =  arrayfun(@(x) x.name, files, 'UniformOutput', 0);
end
if nargin<2
    idx = [1,2];
end
n = length(files);
k=1;
[entities,info] = loadH5Trace(files{k});
m = length(entities(1).data);
t = (0:m-1) * info.dt;
varargout{1}=t;
stim = cell(length(idx),n);
for ii = 1:length(idx)
    varargout{1+ii} = zeros(m,n);
    varargout{1+ii}(:,k) = entities(idx(ii)).data(:);
    stim{ii,k} = entities(idx(ii)).metadata;
end

for k=2:n
    [entities,info] = loadH5Trace(files{k});
    for ii = 1:length(idx)
        varargout{1+ii}(:,k) = entities(idx(ii)).data(:);
        stim{ii,k} = entities(idx(ii)).metadata;
    end
end
varargout{end+1}=stim;
if nargout<length(varargout)
    warning('You specified too less output arrays. Use n_entities+2 as a thumb rule.')
end

