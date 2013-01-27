function [S] = gaussianKernelMovingAverage(X,window)
%
if ~exist('window','var')
    window = 10;
end
% does not work with NaN values so:
X=X(:);
idx = find(isnan(X));
for ii = 1:length(idx)
    if idx(ii)<1
        X(idx(ii)) = 0;
    elseif idx(ii)>length(X)-1
        X(idx(ii)) = 0;
    else
        if ~isnan(X(idx(ii)-1)) & ~isnan(X(idx(ii)+1)) 
            X(idx(ii)) = (X(idx(ii)-1)+X(idx(ii)+1))./2.0;
        elseif ~isnan(X(idx(ii)-1)) & isnan(X(idx(ii)+1)) 
            X(idx(ii)) = (X(idx(ii)-1));
        elseif isnan(X(idx(ii)-1)) & ~isnan(X(idx(ii)+1)) 
            X(idx(ii)) = (X(idx(ii)+1));
        else
            X(idx(ii)) = 0;
        end
    end
end
%disp(find(X==0))
edges = (-4*window:1:4*window);
%Evaluate the Gaussian kernel
kernel = normpdf(edges,0,window);
center = floor(length(edges)/2); 
%Multiply by bin width so the probabilities sum to 1
%kernel = kernel.*resolution; 
%Convolve spike data with the kernel
S = conv(X,kernel);
%Find the index of the kernel center

%Trim out the relevant portion of the spike density
S = S(center:end-center-1); 
if exist('PLOT','var')
    plot(S,'k')
end