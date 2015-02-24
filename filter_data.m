function xf = filter_data(x,fmin,fmax,srate,type)
% FILTER_DATA Wrapper to some basic filters using filtfilt.
%
% xf = FILTER_DATA(x,fmin,fmax,srate,type)
%   Filter types are:
%       ellip      elliptic filter with parameters 2, 0.1, 40
%       butter     4 pole butterworth filter 
%
% Joao Couto 2013

if ~exist('fmin','var');fmin = 300;end
if ~exist('fmax','var');fmax = 3000;end
if ~exist('srate','var');srate = 20000;end
if ~exist('type','var');type = 'ellip';end

switch type
    case 'ellip'
        if isempty(fmin)
            [b,a]=ellip(2,0.1,40,[fmax]*2/(srate),'low');
        elseif isempty(fmax)
            [b,a]=ellip(2,0.1,40,[fmin]*2/(srate),'high');
        else
            [b,a]=ellip(2,0.1,40,[fmin fmax]*2/(srate));
        end
    case 'butter'
        if isempty(fmin)
            [b,a]=butter(4,[fmax]*2/(srate),'low');
        elseif isempty(fmax)
            [b,a]=butter(4,[fmin]*2/(srate),'high');
        else
            [b,a]=butter(4,[fmin fmax]*2/(srate));
        end
        [b,a]=butter(4,[fmin fmax]*2/(srate));
    otherwise
        disp('unknown filter type')
end
xf=filtfilt(b,a,x);

        