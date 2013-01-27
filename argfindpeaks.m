function idx = argfindpeaks(X, threshold, window)
% idx = findPeakIndexes(X, THRESHOLD, WINDOW)
% FINDPEAKINDEXES Extracts the indexes of the peaks in the data with
% threshold crossing.
%       - discards peaks that occur in the first and last window/2 points.
%       - window is the dead time of the detector.
%       - if X is a matrix, argfindpeaks(X) returns a cell array of the
%       indexes (NOT IMPLEMENTED)
%       - Internal variable PLOT can be used for verbose plot (edit
%       findPeakIndexes).

if ~exist('X','var') || isempty(X)
    error('findPeakIndexes: Data must be a vector;')
end
if ~exist('threshold','var'), threshold = 0; end
if ~exist('window','var'), window = 30; end


X = X(:);
warning('off');%'findpeaks:MinPeakHeight');
[~,idx] = findpeaks(X,'MINPEAKHEIGHT',threshold,'MINPEAKDISTANCE',window);
warning('on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code bellow has been working however the above is probably more 
%efficient.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % npeak = 0;
% % xaux = find(X(floor(window/2.0)+1:end-window-1) > threshold) + floor(window/2.0);
% % xaux0 = 0;
% % 
% % idx = nan(1,length(xaux));
% % 
% % for ii=1:length(xaux)
% %     if (xaux(ii) >= xaux0 + window) && ...
% %             (int32(xaux(ii)+floor(window/2))<length(X)) &&...
% %             (int32(xaux(ii)-floor(window/2))>0)
% %         [~, iaux]=max(X(int32(xaux(ii)-floor(window/2)):int32(xaux(ii)+floor(window/2))));        
% %         npeak = npeak + 1;
% %         idx(npeak) = iaux(1) + xaux(ii) - floor(window/2) - 1;
% %         xaux0 = idx(npeak);
% %     end
% % end
% % idx = idx(~isnan(idx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('PLOT','var')
    figure();
    hold on
    plot(X,'ko-')
    if exist('xaux','var')
        plot(xaux,X(xaux),'r.')
    end
    plot(idx,X(idx),'g.')
    axis tight
    plot(xlim,threshold*[1,1],'b--')
    keyboard
end
