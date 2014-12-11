function [cc,mm,ll] = setFigureDefaults()
set(0, ...
    'DefaultFigureColor', [1,1,1], ...              % Figure properties
    'DefaultFigurePaperType', 'A4', ...
    'DefaultFigurePaperUnits', 'centimeters', ...
    'DefaultFigurePaperSize', [18,10], ...
    'DefaultFigurePaperPosition', [0,0,get(0,'DefaultFigurePaperSize')], ...
    'DefaultFigureColor', 'w', ...
    'DefaultFigureInvertHardcopy', 'off', ...
    'DefaultAxesLayer', 'top', ...                  % Axes properties
    'DefaultAxesTickDir', 'out', ...
    'DefaultAxesNextPlot', 'add', ...
    'DefaultAxesUnits', 'normalized', ...
    'DefaultAxesColor', 'none', ...
    'DefaultAxesBox' ,'off',...
    'defaultAxesLinewidth', 0.8, ...
    'DefaultAxesVisible', 'on', ...
    'DefaultLineLineStyle', '-', ...                % Line properties
    'DefaultLineMarker', 'none', ...
    'DefaultLineMarkerSize', 4, ...
    'DefaultLineMarkerFacecolor', [.2, .2, .9], ... % Colors
    'DefaultLineMarkerEdgecolor', [0,0,0], ...
    'DefaultAxesXcolor', [0, 0, 0], ...
    'DefaultAxesYcolor', [0, 0, 0], ...
    'DefaultAxesZcolor', [0, 0, 0], ...
    'DefaultLineColor', [0, 0, 0], ...
    'DefaultAxesFontSize', 8, ...                  % Fonts and text properties
    'DefaultTextFontSize', 8, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextColor', [0, 0, 0], ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultTextVerticalAlignment', 'bottom', ...
    'DefaultTextHorizontalAlignment', 'left',...
    'DefaultFigureVisible','on');

cc = [  21,  69, 152;...
    230,  20,  36;...
    19, 124,  56;...
    239, 104,  28;...
    82,  24, 126;...
    143,  14,  25;...
    163,  31, 129;...
    96,  96,  96;...
    234,  66,  78;...
    196, 186,  87;...
    74, 136, 256;...
    247, 150,  73;...
    140,  79, 155;...
    140,  79, 155;...
    193,  90,  71;...
    204, 103, 165;...
    165, 165, 165;...
    230,  20,  36;...
    19, 124,  56;...
    21,  69, 152;...
    239, 104,  28;...
    82,  24, 126;...
    143,  14,  25;...
    163,  31, 129;...
    96,  96,  96;...
    234,  66,  78;...
    196, 186,  87;...
    74, 136, 256;...
    247, 150,  73;...
    140,  79, 155;...
    140,  79, 155;...
    193,  90,  71;...
    204, 103, 165;...
    165, 165, 165;...
    230,  20,  36;...
    19, 124,  56;...
    21,  69, 152;...
    239, 104,  28;...
    82,  24, 126;...
    143,  14,  25;...
    163,  31, 129;...
    96,  96,  96;...
    234,  66,  78;...
    196, 186,  87;...
    74, 136, 256;...
    247, 150,  73;...
    140,  79, 155;...
    140,  79, 155;...
    193,  90,  71;...
    204, 103, 165;...
    165, 165, 165;...
    230,  20,  36;...
    19, 124,  56;...
    21,  69, 152;...
    239, 104,  28;...
    82,  24, 126;...
    143,  14,  25;...
    163,  31, 129;...
    96,  96,  96;...
    234,  66,  78;...
    196, 186,  87;...
    74, 136, 256;...
    247, 150,  73;...
    140,  79, 155;...
    140,  79, 155;...
    193,  90,  71;...
    204, 103, 165;...
    165, 165, 165]./256;

mm = ['o','x','+','*','s','d','v','^','<','>','p','h','.',...
    '+','*','o','x','^','<','h','.','>','p','s','d','v',...
    'o','x','+','*','s','d','v','^','<','>','p','h','.'];

ll = {'-','--',':','-.','-','--',':','-.','-',...
    '--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};