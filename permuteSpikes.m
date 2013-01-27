function [PERM_SPKTRAIN] = permuteSpikes(SPKTRAIN,OFFSET)
% [PERM_SPKTRAIN] = permuteSpikes(SPKTRAIN)
% SPKTRAIN is an array of doubles.
% Permutes the isis and recalculates the spike times


isis = diff(SPKTRAIN);
isis = isis(randperm(length(isis)));
PERM_SPKTRAIN = ones(size(SPKTRAIN));
randisi = randperm(length(isis));
if ~exist('OFFSET','var')
    PERM_SPKTRAIN(1) = rand(1)*(isis(randisi(1))-isis(randisi(end))/2.0)+...
        SPKTRAIN(1);
else
    PERM_SPKTRAIN(1) = OFFSET;
end

for ii = 1:length(isis)
    PERM_SPKTRAIN(ii+1) = PERM_SPKTRAIN(ii) + isis(ii);
end


