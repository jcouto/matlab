function [info] = mcd2klustersDAT(filename,chunk_size,verbose)
%[info] = mcd2klustersDAT(filename,chunk_size,verbose)
% Converts all analog ('elec') channels from an mcd file to klusters dat format.
% Inputs are the filename (Mandatory), and chunk size (1e6) and verbose (true/false)
if ~exist('filename','var')
    error('Filename not specified.')
end
if ~exist('chunk_size','var')
   chunk_size = 1e6; %floor(N/50);
end
if ~exist('verbose','var')
    verbose = true;
end
log = @(txt)print_if_verbose(txt,verbose);

tstart = tic;
fid = open_mcd_file(filename);

% Build catalogue of entities
% List of EntityIDs needed to retrieve the information and data
[~, nsInfo] = ns_GetFileInfo(fid);
[~, entityInfo] = ns_GetEntityInfo(fid, 1:nsInfo.EntityCount);
%neuralList = find([entityInfo.EntityType] == 4);    
%segmentList = find([entityInfo.EntityType] == 3);
%eventList = find([entityInfo.EntityType] == 1);
channelType = 'elec';
analogList = find([entityInfo.EntityType] == 2);
tmp = cellfun(@(x)x(1:4),{entityInfo(analogList).EntityLabel},'uniformoutput',0);

channelList = analogList(strcmp(tmp,channelType));
channelNumber = cellfun(@(x)str2double(x(14:18)),{entityInfo(channelList).EntityLabel})';
% Create chunks list
N  = ceil(nsInfo.TimeSpan./nsInfo.TimeStampResolution);

nchannels = length(channelNumber);
chunks = 1:chunk_size:N;
if chunks(end) ~= N
    chunks(end+1) = N;
end
log(sprintf('Reading %d samples (%d chunks) from %d analog channels.\\n',...
    N,length(chunks),nchannels));

% Parameters
if ~length(channelNumber)
    error('No analog channels.')
end
entityID = channelList(channelNumber==channelNumber(1));
[~, analogInfo] = ns_GetAnalogInfo(fid, entityID);
info.srate = analogInfo.SampleRate;
info.MCDrange = [analogInfo.MinVal,analogInfo.MaxVal];
info.MCDunits = analogInfo.Units;
info.MCDresolution = analogInfo.Resolution;
range = 10;
amplification = 1000;
info.range = range;
info.amplification = amplification;
info.nBits = 16;
info.nchannels = nchannels;
info.nsamples = N;

% Initialize file and map it to memory
bfilename = strrep(filename,'mcd','dat');
if exist(bfilename,'file')
    log('File already exists, skipping conversion.')
    return 
end
bfid = fopen(bfilename,'w','b');
if N*nchannels*2 > 1e9
    log('File is bigger than 1gb, creating in chuncks\n');
    for i = 1:length(chunks) - 1
        fwrite(bfid,zeros(chunks(i):chunks(i)+diff(chunks(i:i+1))-1,...
            nchannels,'uint16'),'int16');
    end
else
    fwrite(bfid,zeros(N,nchannels,'int16'),'int16');
end
fclose(bfid);

mfile = memmapfile(bfilename,     ...
    'Format', {'int16' [N nchannels] 'data'},  ...
    'Repeat', 1, 'Writable', true);

log('DAT file mapped to memory.\n');
log(['   [',arrayfun(@(z)' ',1:(length(chunks))-1),']\n'])


for i = 1:length(chunks) - 1
    for  j = 1:nchannels
        entityID = channelList(channelNumber==channelNumber(j));
        [~,cont_count, x] = ns_GetAnalogData(fid,entityID,chunks(i),diff(chunks(i:i+1)));
        % Need a conversion to int16
        mfile.Data.data(chunks(i):chunks(i)+cont_count-1,j) = int16((x(:).*amplification*(2^16)/2)./(range./2));
    end
    log([repmat('\b',1,length(chunks)+2),...
        ['[',arrayfun(@(z)'=',1:i),arrayfun(@(z)' ', ...
        1:(length(chunks)-i-1)),']\n']]);
end
ttaken = toc(tstart);
log([repmat('\b',1,length(chunks)-1 + 6),'Conversion toke ',num2str(ttaken),' sec.\n']);

function print_if_verbose(text,verbose)
if verbose
    fprintf(1,text);
end
