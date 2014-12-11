function [G,Ph,freq_edges,param] = transfer_function_corr_analysis(V,I,dt,tL,freq_edges)
% tL - time lag

if ~exist('tL','var');tL = [];end
if ~exist('freq_edges','var');freq_edges = [];end
if isempty(tL)
    tL = 0.5;
end
if isempty(freq_edges)
    freq_edges = 10.^(0:0.01:3.01);
end

time = (0:length(V)-1)*dt;
L = find(time<tL,1,'last');
[Rin,lags] =xcorr(I,L,'biased');

spkidx = argfindpeaks(V,-20);
ts = time(spkidx);
idx = find((ts>tL) & (ts<time(end)-tL));
spkidx = spkidx(idx);
ts = ts(idx);

sta = nan(2*L+1,length(ts));
% spkA = nan(2*L+1,length(ts));
for ii = 1:length(spkidx)
%     spkA(:,ii) = V(spkidx(ii)+[-L:L]);
    sta(:,ii) = I(spkidx(ii)+[-L:L]);
end

tau = lags*dt;

csr = mean(sta,2)';
css = Rin;
Csr = nan(size(freq_edges));
Css = nan(size(freq_edges));

for ii = 1:length(freq_edges)
    f = freq_edges(ii);
    Csr(ii) = sum(csr.*exp(-(tau.^2).*(f.^2)./2).*exp(-1i.*2.*pi.*f.*tau));
    Css(ii) = sum(css.*exp((-(tau.^2).*(f.^2))./2).*exp(-1i.*2.*pi.*f.*tau));
end

G = abs(Csr)./abs(Css);
Ph = atand(imag(Csr)./real(Csr));
[~, tau_delay] = max(csr);
tau_delay = tau(tau_delay);

Ph = Ph-360.*(f*tau_delay);

isi = diff(ts);
param = struct('CV',std(isi)./mean(isi),...
    'frate',mean(isi),...
    'csr',csr,...
    'css',css,...
    'tau',tau);