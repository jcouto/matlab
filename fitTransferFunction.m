function [A,Z,P,Dt,offset,Mopt,Nopt,omega,magnitude,phase,fplot] = ...
    fitTransferFunction(F, magnitude, phase, M, N, weight, coeff, do_plot, robust, guess)
% FITTRANSFERFUNCTION fits a rational transfer function to the modulus and phase of
% frequency response data.
% 
% [A,Z,P,Dt,offset,Mopt,Nopt,omega,magnitude,phase] = ...
%           fitTransferFunction(F, magnitude, phase, M, N, weight, coeff, do_plot, robust, guess)
% 
% The fitted transfer function has the form:
% 
%                M
%              -----
%               | |
%               | |  (jw-Z_k)
%               k=1
%   T(jw) = A ---------------
%                N
%              -----
%               | |
%               | |  (jw-P_k)
%               k=1
% 
% Parameters:
%           F - the frequency values at which the magnitude was measured (in Hertz).
%   magnitude - the magnitude of the response of the system.
%       phase - the phase of the response of the system.
%           M - the degree of the numerator polynomial (i.e., the number of
%               zeros of the transfer function).
%           N - the degree of the denominator polynomial (i.e., the number of
%               poles of the transfer function. N should be greater than M.
%      weight - values used to weigh the error. It must be a L-by-2
%               vector, where L is the number of elements of the frequency
%               vector F. The first and second rows are used for weighing the
%               optimization of the fit of the magnitude and of the phase,
%               respectively. If you don't want to use weights, pass an
%               empty array.
%       coeff - a 2x1 vector of coefficients used to weigh the relative
%               contribution of the magnitude and phase errors in the choice of the
%               best (M,N) combination.
%     do_plot - whether to generate a plot of the results (default is no).
%      robust - whether to prompt the user for the goodness of the fit (default is no).
%       guess - a structure containing three fields, A, Z and P, that were
%               returned by a previous call to fitTransferFunction. In this
%               case, only the optimization of the phase values is
%               performed. Additionally, if guess is present, M and N are
%               automatically taken to be length(A.Z) and length(A.P),
%               respectively.
% 
% NOTE: M and N can either be scalar or vectors: in the latter case, all possible
%       combinations of number of poles and zeros are tested and the function returns
%       the one that gives the smallest (weighted) sum of the errors in the fits of the
%       magnitude and of the phase.
%
% Returns:
%           A - the constant term in the expression of the transfer function.
%           Z - the zeros of the transfer function (in Hertz).
%           P - the poles of the transfer function (in Hertz).
%          Dt - the delay that accounts for the best fit of the phase of
%               the transfer function.
%      offset - the offset in the phase, in the set (-pi,0,pi).
%        Mopt - if M is a vector, Mopt is the number of zeros that leads to the
%               best fit, otherwise Mopt = M.
%        Nopt - if N is a vector, Nopt is the number of poles that leads to the
%               best fit, otherwise Nopt = N.
%       omega - the frequency values at which the magnitude was optimized (in radians/sec).
%   magnitude - the magnitude of the response of the system.
%       phase - the phase of the response of the system (NOTE that the
%               values of phase might be different from those that were
%               passed to the function, if the 'robust' option was set.
% 

% Author: Daniele Linaro - June 2013.

if ~exist('weight','var') || isempty(weight)
    weight = ones(length(F),2);
end
if all(size(weight) == [2,length(F)])
    weight = weight';
end
if any(size(weight) ~= [length(F),2])
    error('weight must be a L-by-2 matrix, with L=length(F).');
end
if ~exist('coeff','var') || isempty(coeff)
    coeff = [1;1];
end
if ~exist('do_plot','var') || isempty(do_plot)
    do_plot = 0;
end
if ~exist('robust','var') || isempty(robust)
    robust = 0;
end
if ~exist('guess','var')
    guess = [];
else
    M = length(guess.Z);
    N = length(guess.P);
    guess.Z = guess.Z*2*pi;
    guess.P = guess.P*2*pi;
end

F = F(:);
magnitude = magnitude(:);
phase = phase(:);
coeff = coeff(:);

% convert to radians/sec
omega = 2*pi*F;

while 1
    
    erropt = 1e12;
    for m=M
        for n=N
            if m < n
                fprintf(1, 'Optimizing with M=%d and N=%d => ', m, n);
                [a,z,p,dt,dphi,err] = doFit(omega,magnitude,phase,m,n,weight,guess);
                if [err.magnitude,err.phase]*coeff < erropt
                    A = a;
                    Z = z;
                    P = p;
                    Dt = dt;
                    offset = dphi;
                    Mopt = m;
                    Nopt = n;
                    erropt = [err.magnitude,err.phase]*coeff;
                else
                    fprintf(1, 'not ');
                end
                fprintf(1, 'new minimum: E = %g*%g + %g*%g = %g.\n', ...
                    coeff(1), err.magnitude, coeff(2), err.phase, [err.magnitude,err.phase]*coeff);
            end
        end
    end

        
        w = logspace(log10(omega(1)),log10(omega(end)),200);
        phi = -offset + zeros(size(w));
        mag = repmat(A,size(w));
        mag = mag * prod(P) / prod(Z);
        for k=1:Mopt
            phi = phi + atan(w/Z(k));
            mag = mag .* sqrt(w.^2 + Z(k)^2);
        end
        for k=1:Nopt
            phi = phi - atan(w/P(k));
            mag = mag ./ sqrt(w.^2 + P(k)^2);
        end
        fplot.w = w;
        fplot.mag = mag;
        fplot.phi = (phi-w*Dt);
        fplot.f = w./(2*pi);
    if do_plot
        figure('visible','on');
        
        subplot(2,1,1);
        semilogx(w,20*log10(mag),'r-','LineWidth',2);
        hold on;
        semilogx(omega,20*log10(magnitude),'ko','MarkerFaceColor','k');
        for i=1:Mopt
            plot(abs(Z(i))+[0,0],ylim,'r--');
        end
        for i=1:Nopt
            plot(abs(P(i))+[0,0],ylim,'g--');
        end
        axis tight;
        axis([omega([1,end])'.*[0.9,1.1],ylim]);
        box off;
        ylabel('20 log_{10} (r_1/r_0) (dB)');
        title(sprintf('M = %d  -  N = %d', Mopt, Nopt));
        
        subplot(2,1,2);
        semilogx(w,phi*180/pi,'r-','LineWidth',2);
        hold on;
        semilogx(omega,(phase+omega*Dt)*180/pi,'ko','MarkerFaceColor','k');
        for i=1:Mopt
            plot(abs(Z(i))+[0,0],ylim,'r--');
        end
        for i=1:Nopt
            plot(abs(P(i))+[0,0],ylim,'g--');
        end
        axis tight;
        axis([omega([1,end])'.*[0.9,1.1],ylim]);
        box off;
        title(sprintf('Dt = %.2f ms', Dt*1e3));
        xlabel('Frequency (radians/s)')
        ylabel('\Phi');
        
        set(gcf,'Color','w','PaperUnits','Inch','PaperPosition',[0,0,6,5]);
    end

    if ~ robust
        break;
    end
    
    str = input('Is the fit ok? y/n [y]: ', 's');
    if isempty(str) || str=='y'
        break;
    end
    
    while 1
        idx = input('Which phi value do you want to change? ');
        if idx < 0
            idx = length(phase)+1+idx;
        end
        if idx == 0 || idx > length(phase)
            error('Index out of bounds.');
        end
        
        dphi = input('By which amount? ');
        phase(idx) = phase(idx) + dphi;
        
        str = input('Change other values of phi? y/n [n] ','s');
        if isempty(str) || str=='n'
            break;
        end
    end
end

Z = Z/(2*pi);
P = P/(2*pi);

end

function [A,Z,P,Dt,offset,err] = doFit(omega, magnitude, phase, M, N, weight, guess)

opt = optimset('TolFun', 1e-9, 'TolX', 1e-6, 'MaxIter', 100000, 'MaxFunEvals', 1000000);

if ~exist('guess','var') || isempty(guess)
    err.magnitude = 1e3;%1e12;
    % optimize the fit of the magnitude response -> get the poles and zeros of
    % the transfer function
    X0 = ones(M+N+1,1);
    X0(1) = 1;
    X0(2:2+M-1) = 0;
    past = [];
    for k=1:70%10
        [X,e] = fminsearch(@(x) costMagnitude(x,M,N,omega,magnitude,weight(:,1)),X0,opt);
        tmp.e = e;
        tmp.X0 = X0;
        past = [past,tmp];
        %if e < err.magnitude
            fprintf(1,'error:%f < %f, randomizing...\n',e,err.magnitude);
            Xopt = X;
            err.magnitude = e;
            X0 = (-1 + 2*rand(size(X0)));%Xopt + (-1 + 2*rand(size(X0)));
        %end
    end
    [~,tmp] = min([past.e]);
    X0 = past(tmp).X0;
    past = [];
    for k=1:10
        [X,e] = fminsearch(@(x) costMagnitude(x,M,N,omega,magnitude,weight(:,1)),X0,opt);
        tmp.e = e;
        tmp.X0 = X0;
        past = [past,tmp];
        %if e < err.magnitude
            fprintf(1,'error:%f < %f, randomizing...\n',e,err.magnitude);
            Xopt = X;
            err.magnitude = e;
            X0 = Xopt + (-1 + 2*rand(size(X0)));
        %end
    end
    
    [~,tmp] = min([past.e]);
    X0 = past(tmp).X0;
    [X,e] = fminsearch(@(x) costMagnitude(x,M,N,omega,magnitude,weight(:,1)),X0,opt);
    err.magnitude = e;
    % constant term, zeros and poles of the transfer function
    A = Xopt(1);
    Z = abs(Xopt(2:2+M-1));
    P = abs(Xopt(2+M:end));
    
else
    
    A = guess.A;
    Z = guess.Z;
    P = guess.P;
    err.magnitude = 0;

end

% theoretically expected phase, with no delay present
theor_phase = zeros(size(omega));
for k=1:M
    theor_phase = theor_phase + atan(omega/Z(k));
end
for k=1:N
    theor_phase = theor_phase - atan(omega/P(k));
end

% optimize the fit of the phase -> get the value of the delay
err.phase = 1e12;
for dphi=[-pi,0,pi]
    [X,e] = fminsearch(@(x) costPhase(x,omega,phase+dphi,theor_phase,weight(:,2)), 0, opt);
    if e < err.phase
        Dt = X;
        offset = dphi;
        err.phase = e;
    end
end

end

function err = costMagnitude(X, M, N, omega, modulus, weight)
A = X(1);
Z = X(2:2+M-1);
P = X(2+M:end);
y = repmat(A,size(omega));
y = y * prod(P) / prod(Z);
for k=1:M
    y = y .* sqrt(omega.^2 + Z(k)^2);
end
for k=1:N
    y = y ./ sqrt(omega.^2 + P(k)^2);
end
err = norm(abs(modulus-y).*weight);
end

function err = costPhase(Dt, omega, meas_phase, theor_phase, weight)
err = norm(abs((meas_phase+omega*Dt) - theor_phase).*weight);
end
