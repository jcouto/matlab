function [Vthr, height, slope, width, ahdp] = extractAPfeatures(X,dt,PLOT)
% [V_THRESHOLD, HEIGHT, SLOPE, WIDTH,AFTER_POTENTIAL] = extractAPfeatures(X, DT, PLOT)
% Acts on rows.
% Units:
%     - X - mV
%     - DT - s
% Extracts:
%     - threshold voltage (V_THRESHOLD)
%     - height of the action potential from the threshold (HEIGHT)
%     - rise velocity of the action potential (SLOPE)
%     - width of the action potential at half height (WIDTH)
%     - after potential from threshold (AFTER_POTENTIAL)
%     - the correction factor given the interpolated time of the action 
% potential. (NOT IMPLEMENTED)
%
if ~exist('dt','var'),dt = 1./15e3;end
if length(dt)>1,dt = 1./(dt(2) - dt(1)); end
if ~exist('PLOT','var'),PLOT = 0;end

t            = (0:size(X,2)-1)*dt*1e3;

INTERP_RATE  = 30e3;
METHOD       = 'spline';
RMSFACTOR    = 0.5;
SMOOTH3dV    = 5; % if >1 the 3rd derivative is smoothed with a moving 
% average of SMOOTH3dV size
ti           = 0 : 1e3/INTERP_RATE : t(end);          % ms
disp(t(end))
dti          = ti(2)-ti(1);

% export variables
N       = size(X,1);
Vthr    = nan(N,1);
slope   = nan(N,1);
height  = nan(N,1);
width   = nan(N,1);
ahdp    = nan(N,1);

%for ii = 1:N % Need to use standard FOR if PLOT is defined.
parfor ii = 1:N
    % cubic spline interpolation
    %disp(['Running ',num2str(ii)])
    
    Xi              = interp1(t,X(ii,:),ti,METHOD);
    [Vmax , iMax]       = max(Xi);                              % !!! mV
    d3Xi            = smooth(diff(Xi,3),SMOOTH3dV,'moving')/((dti*1e3)^3);
    %[~, i3dV]       = max(abs(d3Xi);
    [~,i3dV]        = findpeaks((d3Xi),'minpeakheight',rms(d3Xi)*RMSFACTOR,'minpeakdistance',3);
    
    i3dV = i3dV(find(ti(i3dV)>ti(iMax)-2,1,'first')); %find the first 2ms before the spike
%     iThr            = i3dV + 3; % correct for idx of the 3rd derivative
    iThr            = i3dV;
    Vthr(ii)            = Xi(iThr);                             % !!! mV
    
    slope(ii)           = (Vmax - Vthr(ii))/(ti(iMax)-ti(iThr));    % !!! mV/ms
    height(ii)          = (Vmax - Vthr(ii));                        % !!! mV
    Vwidth              = Vmax-(height(ii)/2.0);
    
    tWidth             = [ti(find(Xi>=Vwidth,1,'first')),... 
                       ti(find(Xi>=Vwidth,1,'last'))];
                   
    width(ii)           = diff(tWidth);                         % ms
    
    d1Xi = diff(Xi,1)/((dti*1e3)^1);
    id1V = find(abs(d1Xi(iMax:end))<rms(d1Xi)/10,1,'first')+iMax;
    ahdp(ii)            = Xi(id1V) - Vthr(ii);                     % !!! mV
    %disp('There is an after hiperpolarization!')

    if PLOT
        ax1 = subplot(2,1,1);
        tWidth = tWidth-ti(iThr);
        tip = ti-ti(iThr);
        hold on
        
        plot(tip(iThr+(-3:3)),Xi(iThr*ones(7)),'r','linewidth',2)
        plot(tip(iMax+(-3:3)),Xi(iMax*ones(7)),'b','linewidth',2)
        plot(tip(id1V(end)+(-3:3)),Xi(id1V(end)*ones(7)),'g-','linewidth',2)
        if N < 2
            text(tip(iThr-5),Vthr(ii),{['slope: ',num2str(slope(ii),'%3.1f'),'mV/ms'],...
                ['height: ',num2str(height(ii),'%2.1f'),'mV/ms'],...
                ['threshold: ',num2str(Vthr(ii),'%2.1f'),'mV'],...
                ['spike width: ',num2str(width(ii),'%2.2f'),'ms'],...
                ['\DeltaAH/DP: ',num2str(ahdp(ii),'%2.2f'),'mV']},...
                'fontname','arial','fontsize',10,...
                'verticalalignment','bottom',...
                'horizontalalignment','right');
        end
        plot([tWidth(1)-0.1,tWidth(2)+0.1],[1,1]*Vwidth,'-','color',[0.5,0.5,0.5],'linewidth',4);
        plot(tip,Xi,'k','linewidth',1)
        axis auto
        set(gca,'box','off','fontname','arial','fontsize',10)
        ax2 = subplot(2,1,2);
        plot(tip(4:end),d3Xi,'r','linewidth',1)
        hold on
        plot(tip(i3dV+3),d3Xi(i3dV),'.k')
        linkaxes([ax1,ax2],'x')
        plot(xlim(),[1,1]*RMSFACTOR*rms(d3Xi),'b--')
        %plot(tip(iMax+(-3:3)),Xi(iMax*ones(7)),'b','linewidth',2)
        set(gca,'box','off','fontname','arial','fontsize',10)
    end
end

   