function V = integrate_eLIF(x,t,I,C,V0,Treset)

if ~exist('Treset','var')
    Treset = 2e-3;
end
if ~exist('V0','var')
    V0 = x(2);
end
f = @(a,X)(1.0/a(1))*(a(2) - X + (a(3) * exp((X - a(4))/a(3))));

V = nan(size(t));
V(1) = V0;

dt = (t(2)-t(1))*1.0e3;
t_last_spk = -10;

for ii = 2:length(I)
    if (t(ii)-t_last_spk) > Treset
        if isnan(V(ii-1))
            % for the fancy plot were the refractory period is blacked
            V(ii-1) = x(2);
        end
        V(ii) = V(ii-1) + dt*(f(x,V(ii-1))+I(ii)/(C));
        if V(ii) > 0
            V(ii) = 20;
            t_last_spk = t(ii);
        end
    end
    
    
end