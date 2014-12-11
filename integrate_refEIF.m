function V = integrate_refEIF(rPar,par,dt,I,C,rFunc,Vreset,tarp)
% integrate_refEIF integrates an exponential integrate and fire
% model with refractoriness.
% 
% V = integrate_refEIF(rPar,dt,I,C,rFunc,Vreset,tarp)
% 
% Parameters:
%   rPar - parameters of the refractory fit, to be used in rFunc (cell array)
%    par - parameters of the EIF model
%     dt - the time step of the simulation, in ms.
%      I - the injected current in pA, in a Nx1 vector.
%      C - the capacitance of the cell, in pF.
%  rFunc - function handles used in the fit (cell array).
%     Vreset - initial condition of the simulation (default, E_m).
%   tarp - absolute refractory period (default, 2 ms).
% 
% Returns:
%      V - a Nx1 vector containing the simulated membrane voltage.

if ~exist('tarp','var')
    tarp = 2;
end
if ~exist('Vreset','var')
    Vreset = [];
end
if isempty(Vreset)
    Vreset = rFunc{2}(par(2),rPar{2},1000);
end
f = @(p,rP,f,tspk,vm)f{1}(p(1),rP{1},tspk)*...
    (f{2}(p(2),rP{2},tspk) - vm + ...
    f{3}(p(3),rP{3},tspk) * ...
    exp((vm - f{4}(p(4),rP{4},tspk))/f{3}(p(3),rP{3},tspk)));

V = Vreset + zeros(size(I));
V(1) = Vreset;
t_since_last_spike = 2*tarp;
tarp_i = round(tarp./dt);
i = 2;
while i <= length(I)
    t_since_last_spike = (t_since_last_spike + dt);
    V(i) = V(i-1) + dt*(f(par,rPar,rFunc,t_since_last_spike,V(i-1)) + I(i-1)/C);
    if V(i) > 0
        V(i) = 40;
        i = i+tarp_i;
        t_since_last_spike = tarp;
    else
        i = i+1;
    end
end
