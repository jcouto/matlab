function indexes = findTransitions(data,thresh)
% FINDTRANSITIONS returns the indexes where a transition occured.
% 
% indexes = findTransitions(data,thresh)
% To be used on sharp transition data like a square waveform.
% Default value for thresh is half of the difference of the max and min
% values of the 1st derivative of the data.
% Joao Couto, June 2012

    diffdata=diff(data);
     if nargin<2 || isempty(thresh)
        thresh=(max(diffdata)-min(diffdata))/2+min(diffdata);
     end
    indexes=find(abs(diffdata)>thresh)+1;

