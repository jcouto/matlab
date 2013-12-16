function V = integrate_EIF(p,dt,I,C,V0,tarp)
% integrate_EIF integrates an exponential integrate and fire
% model.
% 
% V = integrate_EIF(p,dt,I,C,V0,Treset)
% 
% Parameters:
%      p - the vector of parameters of the model. It must be
%          a 1x4 vector, containing the parameters tau_m (ms),
%          E_m (mV), delta_T (ms) and V_T (mV).
%     dt - the time step of the simulation, in ms.
%      I - the injected current in pA, in a Nx1 vector.
%      C - the capacitance of the cell, in pF.
%     V0 - initial condition of the simulation (default, E_m).
%   tarp - absolute refractory period (default, 2 ms).
% 
% Returns:
%      V - a Nx1 vector containing the simulated membrane voltage.
% 

if ~exist('tarp','var')
    tarp = 2;
end
if ~exist('V0','var')
    V0 = [];
end
if isempty(V0);V0=p(2);end

f = @(par,vm) 1/par(1) * (par(2)-vm+par(3)*exp((vm-par(4))/par(3)));
V = p(2) + zeros(size(I));
V(1) = V0;
tarp_i = round(tarp/dt);
i = 2;
while i <= length(I)
    V(i) = V(i-1) + dt*(f(p,V(i-1)) + I(i-1)/C);
    if V(i) > 0
        V(i) = 40;
        i = i+tarp_i;
    else
        i = i+1;
    end
end
