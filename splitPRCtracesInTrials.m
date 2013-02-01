function [trials,t0,trialSpkShapes,baselineCurrent,pertAmp,pertArea] = splitPRCtracesInTrials(T, Vdata,Idata,Pdata,Vthreshold,Pthreshold,nSpkBefore,nSpkAfter)
% This function separates a trace into different trials. It assumes that a
% traditional/direct Phase Response Curve protocol was used.
%
%[trials,t0,spk_shapes,baselineCurrent,pertAmp,pertArea] = splitPRCtracesInTrials(T, Vdata,Idata,Pdata,Vthreshold,Pthreshold)
%     * T can be either the time vector or the dt.
%


if ~exist('Vthreshold','var'), Vthreshold = -10; end
if ~exist('Pthreshold','var'), Pthreshold = 10; end
if length(T)==1, 
    dt=T; T=(0:length(Vdata)-1).*dt;
else
    dt = T(2)-T(1);
end

TPRE = 6;         % (ms)
TPOST = 6;        % (ms)
TDEAD = 5;        % (ms)

Iperturbations = argfindpeaks(Pdata,Pthreshold);
Ispks = argfindpeaks(Vdata,Vthreshold);

Tspks = T(Ispks);
Tperturbations = T(Iperturbations);

N = length(Iperturbations);

trials = cell(1,N);
trialSpkShapes = nan(N,1+TPRE./1000./dt + TPOST./1000.0/dt);
t0 = nan(N, 1);
baselineCurrent = nan(N, 1);
% Amplitude of the perturbation
pertAmp = nan(N, 1);
% Area of the perturbation
pertArea = nan(N, 1);

%DEBUGPLOT = 1;
if exist('DEBUGPLOT','var')
    figure(1),clf
    plot(T,Vdata,'k')
    hold on
    plot(Tspks,Vdata(Ispks),'ko','markerfacecolor','r')
    plot([Tperturbations;Tperturbations],[ones(size(Tperturbations))+10;ones(size(Tperturbations))],'g-')
end
idx = [];
parfor jj = 1:N
    spksToExtract = unique([find(Tspks<=Tperturbations(jj),nSpkBefore,'last'),...
        find(Tspks>Tperturbations(jj),nSpkAfter,'first')]);
    if length(spksToExtract) == (nSpkBefore+nSpkAfter);
        tmpspks = Tspks(spksToExtract) - Tperturbations(jj);
        % In miliseconds
        trials{jj}             = tmpspks*1.0e3;
        % Shape of the spikes
        trialSpkShapes(jj,:)   = nanmean(extractTriggeredTraces(Vdata,Ispks(spksToExtract),...
            TPRE./1000./dt , TPOST./1000./dt));
        % Absolute timestamp of perturbation (MATLAB format)
        t0(jj) = Tperturbations(jj);
        % Current injected (e.g. PID)
        injectedCurrent(jj) = nanmean(extractTriggeredTraces(Idata,Iperturbations(jj),...
            TPRE./1000./dt , TPOST./1000./dt));
        
        pert_traces = extractTriggeredTraces(Pdata,Iperturbations(jj),...
            TPRE./1000./dt , TPOST./1000./dt);
        pertAmp(jj) = nanmax(pert_traces)-nanmin(pert_traces);
        pertArea(jj) = (trapz(pert_traces)*dt); % in miliseconds later this is divided by perturbation amplitude.
        % Mean firing rate (Trial)
        %frate(jj) = 1./nanmean(diff(tmpspks));
    else
        idx = [idx,jj];
    end
end
if length(idx)
    trials(idx) = [];
    t0(idx) = [];
    trialSpkShapes(idx,:) = [];
    baselineCurrent(idx) = [];
    pertAmp(idx) = [];
    pertArea(idx) = [];
end
