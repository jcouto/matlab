function xf = filter_data(x,fmin,fmax,srate,type)
% xf = filter_data(x,fmin,fmax,srate,type)
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

        