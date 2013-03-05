function [phi, dphi, expected_isi, pidx]=extractPRC(spiketimes,option)
% This function extract the raw phase response curves from a spiketimes
% array (each cell is a set of spiketimes (that have zero as the 
% time of the perturbation for each set).
%
% [phi, dphi, expected_isi, idx]=extractPRC(spiketimes,optional)
%
% Optional inputs:
%     * optional is a string that sets the method used to compute the mean
% interspike interval. [unperturbed_mean|unperturbed|mean]
% Outputs:
%     * phi and dphi are the prc points for each order. The minimum order
%     possible for all trials is taken.
%     * idx is the index of the prc that has order zero.
%     * mean_isi is the <isi> considered to calculate dphi
% Notes:
%     * the values are not normalized to the mean isi.
% Joao Couto - February 2013

% Index of perturbed isi. 
idx = cellfun(@(x) find(x>0,1)-1,spiketimes,'UniformOutput',0);
isis = cellfun(@(x) diff(x),spiketimes,'UniformOutput',0);

% Calculate the <isi> that is used to compute dphi default method 'unperturbed_mean'
if ~exist('option','var'), option = 'unperturbed_mean'; end
switch option
    case 'unperturbed_mean'
        % <isi> is the mean of the unperturbed isis.
        expected_isi = cellfun(@(x,y)mean(x(1:y-1)),isis,idx,'UniformOutput',0);
    case 'unperturbed'
        % <isi> is the isi that preceds the perturbation.
        expected_isi = cellfun(@(x,y)mean(x(y-1)),isis,idx,'UniformOutput',0);
    case 'mean'
        % <isi> is the mean of all isis
        expected_isi = cellfun(@(x)mean(x),isis,'UniformOutput',0);
    otherwise
        error('Specified method/option not defined.')
end
% Temporary cell arrays to hold the phase and delta phase points.
phi_ = cellfun(@(x,y)i_calculate_phi(x,y),spiketimes,idx,'UniformOutput',0);
dphi_ = cellfun(@(x,y)i_calculate_dphi(x,y),isis,expected_isi,'UniformOutput',0);

% Find out the minimum number of orders to the left and to the right.
Ml = min(cellfun(@(x,y)length(x(1:y-1)),phi_,idx));
Mr = min(cellfun(@(x,y)length(x(y+1:end)),phi_,idx));
% Prepare outputs
phi = cell2mat(cellfun(@(x,y)x(y-Ml:y+Mr),phi_,idx,'UniformOutput',0))';
dphi = cell2mat(cellfun(@(x,y)x(y-Ml:y+Mr),dphi_,idx,'UniformOutput',0))';
expected_isi = cell2mat(expected_isi);
pidx = Ml+1;


function p = i_calculate_phi(spks, idx)
% Internal function of "extractPRC"
% Computes phi for one trial and all orders.
% Inputs:
%     * spks - spike timestamps
%     * idx - the index of the perturbed isi
%
M = length(spks)-1;
p = nan(M,1);
isi = diff(spks);
for ii = -(idx-1):M-idx
    % To be sure lets do one order at a time.
    if ii == 0
        p(ii+idx) = abs(spks(idx));
    elseif ii < 0
        p(ii+idx) = sum(isi(ii+idx:idx-1))+abs(spks(idx)); 
    elseif ii > 0
        
        p(ii+idx) = abs(spks(idx))-sum(isi(idx+1:idx+ii));
    end
end

function dp = i_calculate_dphi(isi, predicted_isi)
% Internal function of "extractPRC"
% Computes the dphi for one trial and all orders.
% Inputs:
%     * isi - interspike intervals
%     * idx - the index of the perturbed isi
%
dp = predicted_isi-isi(:);
