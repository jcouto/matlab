function out=plotRastergram(varargin)

    p=inputParser;
    p = inputParser;
    
    p.addRequired('spikes') %spike times in secconds
    p.addOptional('counter',0,@isnumeric) %spike times in secconds
    p.addOptional('linewidth',1,@isnumeric) %spike times in secconds

    p.addOptional('marker','o')
    p.addOptional('color',[0,0,0])
    p.addOptional('type','few',@(x)strcmpi(x,'massive')|| strcmpi(x,'few'))
    
    
    
    p.parse(varargin{:})
    spikes=p.Results.spikes;
    type=p.Results.type;
    counter=p.Results.counter;
    marker=p.Results.marker;
    cc=p.Results.color;
    linesize=p.Results.linewidth;
    
    out=gca;
    
    if ~iscell(spikes)
        disp('---> Spike Times must be a cell array...')
        return
    end
    
    if min(size(spikes))>1
        disp('---> Only one dimentional cell arrays can be plotted...')
        return
    end
    hold on
    switch type
        case 'massive' 
            for ii=1:length(spikes)
                plot(spikes{ii},ones(size(spikes{ii})).*counter,marker,'color',cc)
                counter=counter+1;
            end
            
        case 'few' 
            minimum=-0.4;
            maximum=0.4;
            for ii=1:length(spikes)
                for jj=1:length(spikes{ii})
                    line(spikes{ii}(jj).*[1,1],[minimum,maximum]+counter,'color',cc,'linewidth',linesize)
                end
                counter=counter+1;
            end
        otherwise
            disp('---> What is it that you want? Options are (massive) and (few)...')
            return
    end
    %xlabel('Time (s)','fontname','Arial','fontsize',12,'fontweight','bold')
    %ylabel('Trials','fontname','Arial','fontsize',12,'fontweight','bold')