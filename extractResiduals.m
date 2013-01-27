function [xR, R, binnedR] = extractResiduals(x0, y0, xm, ym)
%[xR, R, binnedR] = extractResiduals(x0, y0, xm, ym)
%
% Returns the residuals at each point of x0 using the closest point of the
% model data.
%
%       - x0 and y0 are "measured" data points
%       - xm and ym are the model of the dataset.
%       - if the x0 point is farther than the dt of xm the residual is not
%       calculated.

% remove NaN
inv=(x0~=x0)|(y0~=y0);
x0(inv)=[];
y0(inv)=[];
%
R = nan(length(x0),1);
xR = nan(length(x0),1);

binnedR = cell(length(xm),1);
dt = xm(2)-xm(1);

N = length(x0);
VERBOSE = 0;
for ii = 1:N
    idx = find(xm<=x0(ii),1,'last');
    if isempty(idx)
        idx = 1;
    end
    if (abs(diff([xm(idx),x0(ii)]))<=dt)
        R(ii) = y0(ii)-ym(idx);
        % valid residuals for each x0
        xR(ii) = x0(ii);
        % valid residuals for each point of xm
        binnedR{idx} = [binnedR{idx},R(ii)];
    else
        if VERBOSE
            disp(['Point ',num2str(x0(ii)),' out of reach ',num2str(xm(idx)),'.'])
        end
    end % if
end %for