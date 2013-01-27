function traces = extractTriggeredTrace(x, idx, npre, npost)
% EXTRACTTRIGGEREDTRACE: Extracts the traces given a trigger vector.
%
% NPRE is the number of points to extract before the trigger and NPOST the
% points after.
%
% TRACES that can not be returned within (X,WPRE,WPOST) are returned as
% NAN.
%
% Usage: TRACES = extractTriggeredTrace(X, IDX, NPRE, NPOST)

% handle inputs
if ~exist('x','var') || ~exist('idx','var'),
    error('extractTriggeredTraces: You must specify x and idx.');
end
if ~exist('npre','var'), npre = 40;end
if ~exist('PLOT','var'), PLOT = 0;end
if ~exist('npost','var'), npost = 40;end

% code
if isempty(idx)
    traces = [];
    return
else
    traces = nan(length(idx),npre+npost+1);
end
parfor ii=1:length(idx)
    if ((idx(ii)+npost+2)<=(length(x)-1)) && ((idx(ii)-npre-1)>=1) % check if within boundaries
        traces(ii,:)=x(idx(ii)-npre:idx(ii)+npost);  
    end
end

%%%%%%%%%%%%%%%%%%%%% USELESS %%%%%%%%%%%%%%%%%%%%%
PLOT = 0;
if PLOT 
    fig = figure();axes();
    set(gca, 'ColorOrder', summer(length(idx)), 'NextPlot', 'replacechildren','color','k');
    step = 3;
    plot3( ...
        repmat([0:1:1*size(traces,2)-1]',1,size(traces,1))',...
        repmat([0:size(traces,1)-1]',1,size(traces,2)), ...
        traces,...
        'linewidth',1)
    view(80,30)
    fnchandle = @(x,y)set(gcf,'userdata',0);
    set(fig,'userdata',1,'CloseRequestFcn',fnchandle);
    zlim([min(min(traces))-1,max(max(traces))+1])
    for ii =0:step:length(idx)-step-1
        if get(gcf,'userdata')
            xlim([ii,ii+step])
            pause(0.01)
        end
    end
    if ~get(gcf,'userdata')
        delete(fig);
    else
        fnchandle = @(x,y)delete(gcf);
        set(fig,'userdata',1,'CloseRequestFcn',fnchandle);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%