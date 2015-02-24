function [info] = nev2klustersDAT(filename,channels,chunk_size,force,verbose)
% NEV2KLUSTERSDAT Converts (Ripple) NEV files to the format used by the Klusters Suite.
%
%   [info] = NEV2KLUSTERSDAT(filename,channels,chunk_size,verbose)
% Converts all (or a selected group of) analog ('elec') channels from an nev file to klusters dat format.
% Inputs are the filename (Mandatory), and chunk size (1e6) and verbose (true/false)
% Output filename is the same as the NEV but without extension.
%

range = 20;
amplification = 1000;
if ~exist('force','var')
    force = [];
end
if isempty(force)
    force = true;
end

info = [];


if ~exist('filename','var')
    error('Filename not specified.')
end
if ~exist('channels','var')
    channels = [];
end
if ~exist('chunk_size','var')
   chunk_size = 1e6; %floor(N/50);
end
if ~exist('verbose','var')
    verbose = true;
end
log = @(txt)print_if_verbose(txt,verbose);

tstart = tic;
[res,fid] = ns_OpenFile(filename);
% Build catalogue of entities
% List of EntityIDs needed to retrieve the information and data
channelList = find([fid.Entity.ElectrodeID]<1000 & ...
    strcmp({fid.Entity.EntityType},'Analog') & [fid.Entity.FileType] == 3);
if ~isempty(channels)
    channelList = arrayfun(@(x)find(x == [fid.Entity.ElectrodeID] & ...
        strcmp({fid.Entity.EntityType},'Analog') &...
        [fid.Entity.FileType] == 3 ),channels);
end

N  = ceil(fid.TimeSpan);
nchannels = length(channelList);
chunks = 1:chunk_size:N;
if chunks(end) ~= N
    chunks(end+1) = N;
end
log(sprintf('Reading %d samples (%d chunks) from %d analog channels.\\n',...
    N,length(chunks),nchannels));

% Parameters
if isempty(channelList)
    error('No analog channels.')
end
[~, analogInfo] = ns_GetAnalogInfo(fid, channelList(1));

info.srate = analogInfo.SampleRate;
info.NEVrange = [analogInfo.MinVal,analogInfo.MaxVal];
info.NEVunits = analogInfo.Units;
info.NEVresolution = analogInfo.Resolution;

info.range = range;
info.amplification = amplification;
info.nBits = 16;
info.nchannels = nchannels;
info.nsamples = N;
toint16 = str2func(['@(y) int16(((y./',num2str(amplification),').*(2^16))./',num2str(range),')']);
tofloat = str2func(['@(y)((double((y).*',num2str(range),')./(2^16)).*',num2str(amplification),')']);
info.convertToInt16 = toint16;
info.convertToFloat = tofloat;

% Initialize file and map it to memory
bfilename = strrep(filename,'nev','dat');
if exist(bfilename,'file') & ~force
    log('File already exists; skipping conversion.')
    return
end
bfid = fopen(bfilename,'w','b');
if N*nchannels*2 > 1e9
    log('File is bigger than 1gb, creating in chuncks\n');
    for i = 1:length(chunks) - 1
        fwrite(bfid,zeros(diff(chunks(i:i+1))+1,...
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
        [~,cont_count, x] = ns_GetAnalogData(fid,channelList(j),chunks(i),diff(chunks(i:i+1)));
        % Need a conversion to int16
        try
            mfile.Data.data(chunks(i):chunks(i)+cont_count-1,j) = toint16(x(:));
        catch
            keyboard
        end
    end
    log([repmat('\b',1,length(chunks)+2),...
        ['[',arrayfun(@(z)'=',1:i),arrayfun(@(z)' ', ...
        1:(length(chunks)-i-1)),']\n']]);
end
ttaken = toc(tstart);
log([repmat('\b',1,length(chunks)-1 + 6),'Conversion toke ',num2str(ttaken),' sec.\n']);
ns_CloseFile(fid);
save(strrep(filename,'nev','info.mat'),'-struct','info')

function print_if_verbose(text,verbose)
if verbose
    fprintf(1,text);
end