function [dI, dV, F, ADAPT, absI] = analiseLongStepTrace(STIM, DATA, T, THRESHOLD)
% ANALISELONGSTEPTRACE
% Analises the voltage deflection in trace (DATA) at steady state or the frequency response
% to a pulse of current (STIM). T|SRATE can be the time vector or the sampling
% rate.
% [dI, dV, F, I] = analiseLongStepTrace(STIM, DATA, T | SRATE [15e3])
% If multiple protocols are present the last is selected. Internal
% variables passed to find condition: "FIND_NUMBER"
% The function outputs are:
%  * dI     - delta current applied. [0 pA]
%  * dV     - voltage deflection [nan pA]
%  * F      - steady state firing frequency [nan pA]
%  * ADAPT  - adaptation index. The last ISI divided by the first. [nan]
%  * absI   - total current in the step. [0 pA]
%
STIM        = STIM(:);
DATA        = DATA(:);
FIND_NUMBER = -1;  % which protocol to use (-1 or 'end' to specify the last one).
USE_FULL_TRACE_BEFORE = 1;
if ~exist('THRESHOLD','var'), THRESHOLD = -10;end % threshold for spike detection.

STEADY_STATE_TIME = 0.5;    % fraction of tstep to use on the
% calculations of steady state.

USE_FULL_TRACE_BEFORE = 1;

if ~exist('T','var')
    dt = 1./15e3;
    T  = (0:length(DATA)-1).*dt;
else
    if length(T)>1
        dt = T(2)-T(1);
    else
        dt = 1./T;
        T  = (0:length(DATA)-1).*dt;
    end
end


% find transitions
protTransitions = findTransitions(STIM);
tTransitions    = [0 , T(protTransitions)];

idx = find(diff(tTransitions) > 0.05);

if  isempty(idx)  % seconds

    warning('There are no long steps in this protocol. /n considering spontaneous...')
    
    absI    = mean(STIM);
    dI      = 0;
    dV      = 0;
    F       = nan;
    ADAPT   = nan;
    return
    
end
if (ischar(FIND_NUMBER) && strcmp(FIND_NUMBER,'end') || FIND_NUMBER == -1)
    idx = idx(end);
else
    idx = idx(FIND_NUMBER);
end

tstep = tTransitions(idx:idx+1);

if USE_FULL_TRACE_BEFORE
    tstepBefore = tTransitions(find(tTransitions<tstep(1),1,'last'));
else
    tstepBefore = tTransitions(idx) - STEADY_STATE_TIME*diff(tstep);
end
% Does the cell fire?
% Are there spike between the transitions?


dV      = nan;
F       = nan;
ADAPT   = nan;

absI = mean(STIM(T > (tstep(end)-STEADY_STATE_TIME*diff(tstep)) & ...
    T<(tstep(end))));

dI   = absI - mean(STIM(T > tstepBefore &...
    T<(tstep(1))));

spks = T(argfindpeaks(DATA,THRESHOLD));
%PLOT = 1;


if exist('PLOT','var')
    plot(T,DATA,'color',[0.5,0.5,0.5])
    hold on
    plot(T(T > (tstep(end)-STEADY_STATE_TIME*diff(tstep)) & T<(tstep(end))), ...
        DATA(T > (tstep(end)-STEADY_STATE_TIME*diff(tstep)) & T<(tstep(end))),'r')
    plot(T(T > tstepBefore & T<(tstep(1))), ...
        DATA(T > tstepBefore & T<(tstep(1))),'k')
end


if isempty(spks) || isempty(spks(spks>tstep(1)&spks(tstep(2))))
    % then it is a vi trace.
    VSS = mean(DATA(T > (tstep(end)-STEADY_STATE_TIME*diff(tstep)) & T<(tstep(end))));
    dV = VSS - mean(DATA(T > tstepBefore & T<(tstep(1))));
    if exist('PLOT','var')
       disp(['Voltage after: ',num2str(VSS),' , Voltage before: ',num2str(-dV+VSS)]) 
    end
    
else
    spksSS = spks(spks > (tstep(end)-STEADY_STATE_TIME*diff(tstep)));
    F = 1./mean(diff(spksSS));
    if length(spks)>1
        ADAPT = (spks(end) - spks(1))/spks(1);
    end
end
