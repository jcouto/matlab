function ax = axes_grid(dim,width,height,offset,percent,varargin)
% AXES_GRID Create a homogeneous grid of figure axes
%
% ax = AXES_GRID(DIM,WIDTH,HEIGHT,OFFSET,PERCENT,AXES_PROPERTIES)
%
% Creates an axis grid and returns a matrix with the handles to each axes.
%   - DIM: [NROWS NCOLUMNS]
% Optional inputs:
%   - WIDTH: The width used by the axes grid (default 0.8 of figure)
%   - HEIGHT: The height used by the axes grid (default 0.8 of figure)
%   - OFFSET: The offset where the grid starts (default [0.1,0.1] of the figure)
%   - PERCENT: Percent of the space that is used by the axes (default [0.9,0.9] of the axes dimension)
%   - AXES_PROPERTIES: Additional properties passed to the axes uppon
%   creation
%
% Example 1:
%   ax = axes_grid([10,5],.8,0.5,[],[],'box','on');
%   axes(ax(1,2))   % select axes (1,2)
%   plot(1:40)      % Plot to that axes
%
% Creates a grid of 10 rows by 5 columns using 0.8 of the figure height and
% 0.5 of the width. Passes the box:on option to the axes.
%
% Example 2:
%   fig = figure(1);clf
%   ax = axes_grid([10,10],0.9,0.9,[0.05,0.05],1);
%   for i=1:10
%       for j = 1:10
%           axes(ax(i,j))
%           plot(rand(1,10))
%           axis tight
%       end
%   end
%   set(ax,'color','none','box','on','xtick',[],'ytick',[])
%   set(ax(1:end,1),'ytick',[0:0.5:1])
%   set(ax(end,1:end),'xtick',[1:2:10])
%   linkaxes(ax)
%
% Joao Couto (jpcouto@gmail.com) January 2015


% Validate dimensions
if ~exist('dim','var')
    error('Please specify the dimensions of the axes grid.')
end
if isempty(dim)
    error('Please specify the dimensions of the axes grid [NROW NCOL].')
end
% Validate remaining inputs    
if ~exist('width','var');width = [];end
if ~exist('offset','var');offset = [];end
if ~exist('height','var');height = [];end
if ~exist('percent','var');percent = [];end

if isempty(width),width=0.8;end
if isempty(height),height=0.8;end
if isempty(offset),offset=[0.1,0.1];end
if isempty(percent),percent=[0,0]+0.9;end
if length(percent)<length(dim),percent = percent + 0.*dim;end

% Initialize variables and create axes.
h = height./dim(1);
w = width./dim(2);
o = fliplr(offset);
o(1) = o(1) + height;

ax = nan(dim);
for i = 1 : dim(1)
    o = o - [h, 0];
    o(2) = offset(1);
    for j = 1 : dim(2)
        ax(i,j) = axes('position',[o(2), o(1), w*percent(2), h*percent(1)],...
            varargin{:});
        o = o + [0, w];
    end
end