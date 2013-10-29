function [x, f, resnorm, r,o] = fit_eLIF_to_dIV(X, Y, C, x0)
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
%label = {'Time constant (\tau_m)','Resting potential (E_m)', 'Spike width (/delta_T)','Spike threshold (V_t)'}
%%

f = @(a,X)(1.0/a(1))*(a(2) - X + (a(3) * exp((X - a(4))/a(3))));
% Initial values
if ~(exist('x0','var')),x0 = [20,-60,1,-40];end
% Lower bounds
lb = [0,-90,0,-90];
% Upper bounds
ub = [100,-30,15,-30];
options = optimset('TolFun',1e-6,'maxiter',100,'display','off');
[X,idx] = sort(X);
Y = Y(idx);
% figure(9),clf
%  plot(X,-Y/C,'ko')
[x,resnorm,r,~,o] = lsqcurvefit(f,x0,X,-Y/C,lb,ub,options);
%  hold on;plot(X,f(x,X),'r')
%  pause


