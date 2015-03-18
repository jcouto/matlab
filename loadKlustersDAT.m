function out = loadKlustersDAT(filename,nchannels,mode)
% LOADKLUSTERSDAT load klusters DAT file.
% out = LOADKLUSTERSDAT(filename,nchannels,mode)
file = dir(filename);
tmp = int16(1);
s=whos('tmp');
tmps = s.bytes;

nsamples = file.bytes/(tmps*nchannels)
 %for i = 1:nchannels
    dataformat =  {'int16',[nchannels,nsamples],sprintf('lfp%03d',0)};
%end

file = memmapfile(filename,...
    'Format',dataformat)
