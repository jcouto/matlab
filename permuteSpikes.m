function [PERM_SPKTRAIN] = permuteSpikes(spiketrain,N,offset)
% [PERM_SPKTRAIN] = permuteSpikes(SPKTRAIN,OFFSET,N)
% SPKTRAIN is an array of doubles.
% Permutes the isis and recalculates the spike times
SPKTRAIN = spiketrain;
for ii = 1:N
    isis = diff(SPKTRAIN);
    isis = isis(randperm(length(isis)));
    PERM_SPKTRAIN = ones(size(SPKTRAIN));
    randisi = randperm(length(isis));
    if ~exist('OFFSET','var')
        %         PERM_SPKTRAIN(1) = rand(1)*(isis(randisi(1))-isis(randisi(end))/2.0)+...
        %             SPKTRAIN(1);
        PERM_SPKTRAIN(1) = SPKTRAIN(1);%+(rand(1)*0.03*(mean(isis)))-(0.03*(mean(isis))/2);
    else
        PERM_SPKTRAIN(1) = offset;
    end
    
    for ii = 1:length(isis)
        PERM_SPKTRAIN(ii+1) = PERM_SPKTRAIN(ii) + isis(ii);
    end
    SPKTRAIN = PERM_SPKTRAIN;
end
