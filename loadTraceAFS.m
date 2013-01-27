function [out] = loadTraceAFS(varargin)
    % [out]=loadTrace(filename,start=0,stop=0) loads a trace saved with the AFS 
    % input filename can be a string, set of strings or an output of "dir".
    p=inputParser;
    p = inputParser;
    p.addRequired('fname') %filename
    p.addOptional('start',0); %spikes detection threshold
    p.addOptional('stop',0);
    p.parse(varargin{:});

    fname=p.Results.fname;
    Tstart=p.Results.start;
    Tstop=p.Results.stop;
    out=[];
	if isstr(fname) %single file
		out=loadSingleFile(fname,Tstart,Tstop);
	elseif iscell(fname)
		for pp=1:length(fname)
			disp(['--> Reading file: ',fname{pp}])
			tmpout=loadSingleFile(fname{pp},Tstart,Tstop);
			out=[out,tmpout];
		end
	elseif isstruct(fname) %an entire directory
		%then it should be an output of "dir"
		fname={fname.name};
		for pp=1:length(fname)
			disp(['--> Reading file: ',fname{pp}])
			tmpout=loadSingleFile(fname{pp},Tstart,Tstop);
			out=[out,tmpout];
		end
	
	else
		 disp('--> Couldn t really do much read...')
	end
end	


function [out]=loadSingleFile(fname,Tstart,Tstop)
	% is the file there??
	if ~exist(fname, 'file')
        	disp(['---) Warning --> error: file ',fname ,' not found!']);
		out=[];
		%out.ch1=[];
        	%out.srate=-1;
    else
    
    fp = fopen(fname, 'r');        
    srate = fread(fp, 1, 'double');
    N     = fread(fp, 1, 'ulong');
    M     = fread(fp, 1, 'ulong');
    %N     = fread(fp, 1, 'uint64'); % number of simultaneous waves (channels)
    %M     = fread(fp, 1, 'uint64'); % number of samples per wave (channel)
    
    if ~(Tstart==0)
        Tstart=int32(Tstart*srate)+1;
    end
    if Tstop==0
        Tstop=M;
    else
        Tstop=int32(Tstop*srate)+1;
    end
    if Tstop>M
        Tstop=M;
    end
    
    recorded  = zeros((Tstop)-Tstart,N);
    zeroPosition=ftell(fp);
    
    for ii=1:N
        fseek(fp,zeroPosition+(Tstart)*8,'bof'); % move to the correct position (8 is the number of bytes in a double)
        recorded(:,ii)=fread(fp,(Tstop)-Tstart,'double');
        zeroPosition=zeroPosition+(int32(M)*8);
    end
    fclose(fp);
    %export stuff
    out.srate=srate;
    for ii=1:size(recorded,2)
        out.(['ch',num2str(ii)])=recorded(:,ii);
    end
    out.filename=fname;
    end
end
