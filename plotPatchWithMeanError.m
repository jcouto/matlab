function p = plotPatchWithMeanError(X, MEAN, ERROR, FACECOLOR, ALPHA, EDGECOLOR)
% p = plotPatchWithMeanError(X, MEAN, ERROR)
%
% This function is a wrapper around "patch".
% Plots the shadow of the MEAN plus minus the ERROR;

%%
if ~exist('EDGECOLOR','var')
    EDGECOLOR = 'none';
end
if ~exist('FACECOLOR','var')
    FACECOLOR = [.5,.5,.5];
end

if ~exist('ALPHA','var')
    ALPHA = .5;
end

X=X(:);
MEAN = MEAN(:);
ERROR = ERROR(:);

p = fill([X;flipud(X)],[MEAN-ERROR;flipud(MEAN+ERROR)],1,...
    'edgecolor',EDGECOLOR, 'facecolor',FACECOLOR, 'facealpha',ALPHA);
% to ensure that the objects appear in the order they are plot.
uistack(p,'bottom') 
% sets the axis ticks on top of the patch.
if exist('FIXAXES','var')
    set(gca,'DataAspectRatio',[1 1 1],'Layer','top')
end
