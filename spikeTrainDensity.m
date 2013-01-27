function [s, t] = spikeTrainDensity(ts, resolution, sigma)
% [s, t] = spikeTrainDensity(ts, resolution, sigma)
% Returns the spike density function from a spike raster.
%   * ts: spike train (s)
%   * resolution|dt: spike train resolution - binning size (ms)
%   * sigma: confidence interval of the spike times (ms)
% Adapted from Matlab for Neuroscientists.
%

if ~exist('sigma','var')
    sigma = 15; %Standard deviation of the kernel = 15 ms
end
% Bin spike train.
if ~exist('resolution','var')
    resolution = 1; % 1ms bins
end
% Since we want to work in seconds
resolution = resolution/1e3;
sigma = sigma/1e3;
% Use hist to bin
EDGES = (min(ts)-100*resolution:resolution:max(ts)+100*resolution);
N = histc(ts, EDGES);
%Time ranges form -3*st. dev. to 3*st. dev.
edges = (-3*sigma:resolution:3*sigma);
%Evaluate the Gaussian kernel
kernel = normpdf(edges,0,sigma);
%Multiply by bin width so the probabilities sum to 1
kernel = kernel.*resolution; 
%Convolve spike data with the kernel
s = conv(N,kernel);
%Find the index of the kernel center
center = ceil(length(edges)/2); 
%Trim out the relevant portion of the spike density
s = s(center:end-center-1); 
t = (0:length(s)-1)*resolution - abs(min(ts))-100*resolution;
if exist('PLOT','var')
    plot(t,s,'k')
    plotRastergram({ts},'color','r')
end