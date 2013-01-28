function [dI, dV, F, ADAPT, absI] = analiseLongStepTrace(STIM, DATA, ...
                                                  T, THRESHOLD, FRACTION_OF_TSTEP)
% ANALISELONGSTEPTRACE
% Analises the voltage deflection in trace (DATA) at steady state or the frequency response
% to a pulse of current (STIM). T|DT can be the time vector or the
% sampling time
%
% [dI, dV, F, I] = analiseLongStepTrace(STIM, DATA, T | DT [1./30kHz])
%
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

% Which protocol to use (-1 or 'end' to specify the last one).
FIND_NUMBER = -1;  

% Threshold for spike detection.
if ~exist('THRESHOLD','var'), THRESHOLD = -20;end 

% Fraction of tstep to use on the calculations of steady state.
if ~exist('FRACTION_OF_TSTEP','var')
    FRACTION_OF_TSTEP = 0.5;
end

USE_FULL_TRACE_BEFORE = 0;

if ~exist('T','var')
    dt = 1./30e3;
    T  = (0:length(DATA)-1).*dt;
else
    if length(T)>1
        dt = T(2)-T(1);
    else
        dt = T;
        T  = (0:length(DATA)-1).*dt;
    end
end

% Find transitions
protTransitions = findTransitions(STIM);
tTransitions    = [0 , T(protTransitions)];

idx = find(diff(tTransitions) > 0.05);


if  isempty(idx)  
    % Case of the spontaneous protocol.
    warning('There are no long steps in this protocol. Considering spontaneous...')
    
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
    tstepBefore = tTransitions(idx) - FRACTION_OF_TSTEP*diff(tstep);
end
% Does the cell fire?
% Are there spike between the transitions?

dV      = nan;
F       = nan;
ADAPT   = nan;

absI = mean(STIM(T > (tstep(end)-FRACTION_OF_TSTEP*diff(tstep)) & ...
    T<(tstep(end))));

dI   = absI - mean(STIM(T > tstepBefore &...
    T<(tstep(1))));

spks = T(argfindpeaks(DATA,THRESHOLD));

%PLOT = 1;
if exist('PLOT','var')
    plot(T,DATA,'color',[0.5,0.5,0.5])
    hold on
    plot(T(T > (tstep(end)-FRACTION_OF_TSTEP*diff(tstep)) & T<(tstep(end))), ...
        DATA(T > (tstep(end)-FRACTION_OF_TSTEP*diff(tstep)) & T<(tstep(end))),'r')
    plot(T(T > tstepBefore & T<(tstep(1))), ...
        DATA(T > tstepBefore & T<(tstep(1))),'k')
    pause
end
if isempty(spks) || isempty(spks(spks>tstep(1)&spks<tstep(2)))
    % then it is a vi trace.
    VSS = mean(DATA(T > (tstep(end)-FRACTION_OF_TSTEP*diff(tstep)) & T<(tstep(end))));
    dV = VSS - mean(DATA(T > tstepBefore & T<(tstep(1))));
    if exist('PLOT','var')
       disp(['Voltage after: ',num2str(VSS),' , Voltage before: ',num2str(-dV+VSS)]) 
    end
    
else
    spksSS = spks(spks > (tstep(end)-FRACTION_OF_TSTEP*diff(tstep)));
    F = 1./mean(diff(spksSS));
    if length(spks)>1
        ADAPT = (spks(end) - spks(1))/spks(1);
    end
end
