function [x, f, resnorm, r,o] = fit_eLIF_to_dIV(X, Y, C)
% Fits a dynamic IV curve to an exponential integrate and fire model
% neuron.
% It might be necessary to select the values of Y.
% Takes the voltage (X), mean transmembrane current (Y) and the Capacitance (C) as
% parameters
% Returns:
%   - the result of the fit (parameters)
%   - resnorm the squared residuals
%   - f the function used to fit
% The result is organized in the following way: [tau_m, E_m, delta_T, V_T].
% For example:
%idx = (dI_V>-100 & dI_V<=-40 & ~isnan(dI_mu)& -dI_mu/C <= 25 );%
%[x] = fit_dIV_to_eLIF(dI_V(idx), dI_mu(idx),C)
%%

f = @(a,X)(1.0/a(1))*(a(2) - X + (a(3) * exp((X - a(4))/a(3))));
% Initial values
x0 = [20,-60,1,-40];
% Lower bounds
lb = [0,-80,0,-80];
% Upper bounds
ub = [100,-20,15,-20];
options = optimset('TolFun',1e-6,'maxiter',1000);

[x,resnorm,r,~,o] = lsqcurvefit(f,x0,X,-Y/C,lb,ub,options);

