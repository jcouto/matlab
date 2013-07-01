function [A,Z,P,Dt,Mopt,Nopt,plotout] = fitTransferFunction(F, magnitude, phase, M, N, weight, doPlot)
% FITTRANSFERFUNCTION fits a rational transfer function to the modulus of
% frequency response data.
%
% [A,Z,P,Dt,Mopt,Nopt] = fitTransferFunction(F, magnitude, phase, M, N, weight, doPlot)
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
%      weight - values used to weigh the error. It must be a 2-by-L
%               vector, where L is the number of elements of the frequency
%               vector F. The first and second rows are used for weighing the
%               optimization of the fit of the magnitude and of the phase,
%               respectively. If you don't want to use weights, pass an
%               empty array.
%      doPlot - whether to generate a plot of the results (default is no).
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
%        Mopt - if M is a vector, Mopt is the number of zeros that leads to the
%               best fit, otherwise Mopt = M.
%        Nopt - if N is a vector, Nopt is the number of poles that leads to the
%               best fit, otherwise Nopt = N.
%

% Author: Daniele Linaro - June 2013.

if ~exist('weight','var') || isempty(weight)
    weight = ones(2,length(F));
end
if any(size(weight) ~= [length(F),2])
    error('weight must be a 2-by-L matrix, with L=length(F).');
end
if ~exist('doPlot','var')
    doPlot = 0;
end

% convert to radians/sec
omega = 2*pi*F;

erropt = 1e12;
c = [1;0.1];
for m=M
    for n=N
        if m < n
            fprintf(1, 'Optimizing with M=%d and N=%d => ', m, n);
            [a,z,p,dt,err] = doFit(omega,magnitude,phase,m,n,weight);
            if [err.magnitude,err.phase]*c < erropt
                A = a;
                Z = z;
                P = p;
                Dt = dt;
                Mopt = m;
                Nopt = n;
                erropt = [err.magnitude,err.phase]*c;
            else
                fprintf(1, 'not ');
            end
            fprintf(1, 'new minimum: E = %g*%g + %g*%g = %g.\n', ...
                c(1), err.magnitude, c(2), err.phase, [err.magnitude,err.phase]*c);
        end
    end
end
w = logspace(log10(omega(1)),log10(omega(end)),200);
phi = zeros(size(w));
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
plotout.mag = mag;
plotout.phi = phi;
plotout.w = w;
plotout.f = w/(2*pi);

Z = Z/(2*pi);
P = P/(2*pi);

if doPlot
    figure;
    
    
    
    subplot(2,1,1);
    semilogx(w,10*log10(mag),'r-','LineWidth',2);
    hold on;
    semilogx(omega,10*log10(magnitude),'ko','MarkerFaceColor','k');
    for i=1:Mopt
        plot(abs(Z(i))+[0,0],ylim,'r--');
    end
    for i=1:Nopt
        plot(abs(P(i))+[0,0],ylim,'g--');
    end
    axis tight;
    axis([omega([1,end])'.*[0.9,1.1],ylim]);
    box off;
    ylabel('10 log_{10} (r_1/r_0) (dB)');
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

Z = Z/(2*pi);
P = P/(2*pi);

end

function [A,Z,P,Dt,err] = doFit(omega, magnitude, phase, M, N, weight)

err.magnitude = 1e12;
% optimize the fit of the magnitude response -> get the poles and zeros of
% the transfer function
opt = optimset('TolFun', 1e-9, 'TolX', 1e-6, ...
    'MaxIter',10000000,'MaxFunEvals',10000000);
X0 = ones(M+N+1,1);
X0(1) = 1;
X0(2:2+M-1) = 0;
for k=1:10
    [X,e] = fminsearch(@(x) costMagnitude(x,M,N,omega,magnitude,weight(:,1)), X0, opt);
    if e < err.magnitude
        Xopt = X;
        err.magnitude = e;
        X0 = Xopt + (-1 + 2*rand(size(X0)));
    end
end

% constant term, zeros and poles of the transfer function
A = Xopt(1);
Z = abs(Xopt(2:2+M-1));
P = abs(Xopt(2+M:end));

% theoretically expected phase, with no delay present
theor_phase = zeros(size(omega));
for k=1:M
    theor_phase = theor_phase + atan(omega/Z(k));
end
for k=1:N
    theor_phase = theor_phase - atan(omega/P(k));
end

% optimize the fit of the phase -> get the value of the delay
[Dt,err.phase] = fminsearch(@(x) costPhase(x,omega,phase,theor_phase,weight(:,2)), 0, opt);

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
