function [N, L, y, nstim, ntrace, isis] = process_noisy_prc_traces(stimset, vset, dt, meanI, stdI, L, threshold)

%   * N: total number of ISIs,
%   * L: length of each stimuli,
%   * nstim: stretched stimuli for the ISIs in an (NxL) matrix,
%   * y: change in instantaneous firing rate,
%   * ntrace: stretched voltage for the ISIs,
%   * isis: ISI sizes in ms,
%   * ndt: mean ISI size divided by L.
% Adapted from the original code by: Sungho Hong, CNS Unit, Okinawa Inst of Sci Tech
if ~exist('threshold','var'),threshold = -20;end


Nseg = length(vset);
isis = [];

for i=1:Nseg
    v = vset{i};
    spiketime = argfindpeaks(v, threshold);
    isis = [isis; diff(spiketime)];
end

y = (mean(isis) - isis)./isis;
N = length(y);
% length of the rescaled stimuli
% 
if ~exist('L','var') || isnan(L)
    L = round(mean(isis));
end 
if mod(L,2)==0  % L should be always odd.
    L = L+1;
end

D = L;
k = 1;
nstim = zeros(N,L);
ntrace = zeros(N,L);
tt = linspace(0,1,L);
for i=1:Nseg
    stim = (stimset{i} - mean(stimset{i}))/std(stimset{i});
    v = vset{i};
    spiketime = argfindpeaks(v, threshold);
    nisi = length(spiketime)-1;
    for i=1:nisi
        scut = stim(spiketime(i):spiketime(i+1));
        vcut =    v(spiketime(i):spiketime(i+1));
        torig = linspace(0,1,numel(scut));
        nstim(k,:)  = stretch(scut, D); % resample(scut, D, numel(scut)); %interp1(torig,scut,tt, 'pchip');
        ntrace(k,:) = stretch(vcut, D); % resample(scut, D, numel(scut)); %interp1(torig,vcut,tt, 'pchip');
        k = k+1;
    end
end

isis = isis*dt*1.0e3;
dt = mean(isis)/L;

function r = stretch(s_, N_)
% new_signal = stretch(signal, length_of_new_signal) extra/interpoles a row
% vector that preserves the Fourier coefficients. Maybe MATLAB resample
% function does it right, but I want to make sure. 
%
% Inputs
% ======
% signal: original signal
% length_of_new_signal: length of new signal that should be an _odd_
% number since the number of Fourier basis vectors created by dftBasis will
% be alwayas an odd number.
%
% Written by: Sungho Hong, CNS Unit, Okinawa Inst of Sci Tech
% Email: shhong@oist.jp
%

if size(s_,2)==1
    is_s_row = false;
    s = s_;
else
    is_s_row = true;
    s = reshape(s_,[numel(s_) 1]);
end

if mod(N_,2)==0
    error('Length of the new signal should be an odd number for some reasons. Sorry.')
end

L = numel(s);
x = fft(s);
LL = floor(L/2);
N = floor(N_/2);

if (N<=LL-1)
    x = [x(1); x(2:(N+1)); x(end-(N-1):end)];
else
    x = [x(1); x(2:LL); zeros(2*(N-LL)+2,1); x(end-(LL-2):end)];
end

r = ifft(x)*sqrt(N_/L);

if is_s_row
    r = r';
end
