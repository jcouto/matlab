function Vc = AECoffline(Vr,I,Ke)
% Vc = AECoffline(Vr,I,Ke)
% 
% Parameters:
%    Vr - recorded membrane voltage (in mV)
%     I - injected current (in pA)
%    Ke - electrode kernel (in Ohm)
% 
% Returns:
%    Vc - compensated membrane voltage (in mV)
% 

% Daniele Linaro - December 2011

if size(Ke,1) > 1
    Ke = Ke';
end
[r,c] = size(Vr);
if r > c
    Vr = Vr';
end
[r,c] = size(I);
if r > c
    I = I';
end

N = size(Vr,1);
Ue = zeros(size(Vr));
for k=1:N
    tmp = conv(Ke,I(k,:)*1e-12)*1e3;
    Ue(k,:) = tmp(1:end-length(Ke)+1);
end
Vc = Vr - Ue;
