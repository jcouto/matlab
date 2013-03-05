function out = fitExp(varargin)
% FITEXP fits exponentials to data
% 
% Usage:
%       out = fitExp(V,srate,nexp)
% 
% Parameters:
%      V - voltage trace
%  srate - sampling rate
%   nexp - number of exponentials to use in the fit
% 
% Returns:
%    out - structure with F (fit object) and G (goodness-of-fit).

p = inputParser;
p.addRequired('tofit',@isnumeric)
p.addRequired('srate',@isnumeric)
p.addRequired('nexp',@isnumeric)
p.parse(varargin{:});

tofit = p.Results.tofit;
tofit = tofit(:);
nexp = p.Results.nexp;
srate = p.Results.srate;
 

T = (0:length(tofit)-1) ./ srate; %in seconds
T = T(:);

if nexp == 1
     s = fitoptions('Method','NonlinearLeastSquares',...
         'Lower',[-50, 1e-9, 0],...
         'Upper',[50, 1, 15],...
         'Startpoint',[-10, 0.01, 10]);
     ff = fittype('a0 * (1-exp( -x / b0 )) + c',...
         'options',s);
     [F,G] = fit(T,tofit-(min(tofit)),ff);
     %ci = confint(F,0.95);
     out.F = F;%.tau=(ci(:,2));
     out.G = G;%rmse=G.rmse;
 elseif nexp == 2
     %fit options
     s = fitoptions('Method','NonlinearLeastSquares',...
         'Lower',[-50,-50, 1e-9,1e-9, 0],...
         'Upper',[50,50, 1,1, 15],...
         'Startpoint',[5,15, 0.00001,0.005, 10]);
         %'Startpoint',[-10,-10, 0.00001,0.005, 10]);
     ff = fittype('a0 * (1-exp( -x / b0 )) + a1 * (1-exp( -x / b1 )) + c',...
         'options',s);
     %ff = fittype('a0 * exp( -x / b0 ) + a1 * exp( -x / b1 ) + c','options',s);
     
     [F,G] = fit(T,tofit-(min(tofit)),ff);
     %ci = confint(F,0.95);
     out.F = F;
     out.G = G;
 elseif nexp==3
     %fit options
     s = fitoptions('Method','NonlinearLeastSquares',...
         'Lower',[-50,-50,-50, 1e-9,1e-9,1e-9, 0],...
         'Upper',[50,50,50, 1,1,1, 15],...
         'Startpoint',[-10,-10,-10, 0.00001,0.01,0.05, 10]);
     
     ff = fittype('a0 * (1-exp( -x / b0 )) + a1 * (1-exp( -x / b1 )) +  a2 * (1-exp( -x / b2 ))+ c',...
         'options',s);
     
     [F,G] = fit(T,tofit-(min(tofit)),ff);
     ci = confint(F,0.95);
     out.F = F;%tau=(ci(:,4:6));
%      for ii=1:size(out.tau,1)
%          out.tau(ii,:)=sort(out.tau(ii,:));
%      end
     out.G = G;%rmse=G.rmse;
 else
     disp('+++> Set the seccond input to the number of exponentials (1 or 2 in this case...).')
     out.F = [];
     out.G = [];
     return;
end
 
    % b1 and b2 are the values that matter for the resistence of the
    % pipette and the membrane. (the smallest value is for the pipette).
%     
%     if plotflag==1 & ~isempty(F)
%         cla
%         hold on
%         plot(T,waveforms,'ko')
%         %plot(T,mean(waveforms')+std(waveforms'),'color',[1,0.1,0.1])
%         %plot(T,mean(waveforms')-std(waveforms'),'color',[1,0.1,0.1])
%         plot(T,mean(waveforms'),'r','linewidth',2)
%         val=F(T)+min(mean(waveforms'));
%         plot(T,val,'g','linewidth',2);
%         names=[];
%         for ii=1:size(out.tau,1)
%             names(ii,:)=['\tau',num2str(ii),' - '];
%         end
%         text(T(round(length(T)/3)),min(mean(waveforms')),[names, num2str(mean(out.tau)')],'FontWeight','bold')
%         xlabel('time (sec)')
%         ylabel('Voltage (mV)')
%         
%     end