function [freq,psdx,xdft] = compute_power_spectral_density(x,time,PLOT)
% [freq,psdx] = compute_power_spectral_density(x,time,PLOT)
% Computes the Power spectral density
% x is the signal
% time can be only dt.

if ~exist('time','var')
    dt = 1./20e3;
else
    if length(time)>1
        dt = time(2)-time(1);
    else
        dt = time;
    end
end

Fs = 1./dt;

N = length(x);
xdft = fft(x);
if mod(N,2)
    N = N-1;
end

xdft = xdft(2:end-1);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;


if exist('PLOT','var')
    plot(freq,10*log10(psdx),'k'); grid on;hold on;
%     plot(freq,smooth(10*log10(psdx),100),'r');
    axis tight
    %set(gca,'yscale','log')
    axis tight
    xlabel('Frequency (Hz)'); ylabel('Power density (dB)');
    
end