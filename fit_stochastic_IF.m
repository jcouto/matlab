function fit_stochastic_IF(t,V,I)
% This function fits a stochastic IF model to experimental data.
% Adapted from Mensi et al. 2012
%

% convert time to ms and current to nA
t = t*1e3; 
I = I*1e-3;

dt = diff(t(1:2));
%   1. Extract Spike Shape, T_refr and E_reset from V
[spks,spk_wave,tspk_wave, spkidx] = extract_spikes( V, [], t*1e-3, 1, 10, 3);

spks = spkidx*1e3-0.7;
tspk_wave = tspk_wave+abs(tspk_wave(1));
% Set spiketimes to 0.7ms before the actual spike peak (isn't it better to use spike_onsets???)
% figure('visible','on')
% plot(t,V)
% plot(spks-0.0007,V(int32(spkidx-(0.0007/diff(t(1:2))))),'ko')
m_spk_wave = mean(spk_wave);
E_reset = min(m_spk_wave);
temp = find(m_spk_wave==E_reset);                % t_refr is the time after maximum of the AP
t_refr = tspk_wave(temp(end)) - 0.7 ;

size_eta = 2000;                                    % total size of eta in ms
nbr_bink = 40;                                      % nbr bin in eta
temp_bin = @(x,tempdt,bin0) exp(x/tempdt) + bin0;   % function to compute the log-distribution of the kernel bins
tempdt = 0.2;                                       % sharpness of the temp_bin function
x = 0:1/nbr_bink:1-(1/nbr_bink);                    % x-axis of the log-distribution
bink_size0 = 5/dt;                                  % size of the first bin: 5 ms
bin0 = bink_size0*dt;                               % generate the bin distribution and normalize
y = temp_bin(x,tempdt,bin0);
normalization_factor = sum(y); normalization_factor = normalization_factor/size_eta;
bink_size = round(round(y/normalization_factor)/dt);
x=nan(nbr_bink,1);
for j=1:nbr_bink
	if(j==1)
        x(j) = round(bink_size(j)/2.);
    else
        x(j) = round(sum(bink_size(1:j-1))+(bink_size(j)/2.));
	end
end
temp_x = zeros(nbr_bink,sum(bink_size(1:end)));
for j=1:nbr_bink
	temp_x(j,x(j)-round(bink_size(j)/2.)+1:x(j)+round(bink_size(j)/2.)) = ones(1,bink_size(j));
end
% Linear Regression
delay = 1;                                   % parameters used to set the exact timing of the spikes
spike = spks-(delay/dt);             
V = V';
I = I';
diff_v =  [diff(V);0]/dt;                       % compute voltage derivative        
diff_v_int = Extract_interval(diff_v,spike,(delay+t_refr)/dt);  % remove spikes from derivative
M = build_M_matrix(V,spike,delay+t_refr,nbr_bink,bink_size,1000./dt);  % build the matrix that contains the basis function
X = [V ones(length(V),1) I -1*M']';             % build matrix of the regressors
XX = zeros(size(X,1),length(diff_v_int));            % remove spikes from the matrix
for j=1:size(X,1)
	XX(j,:) = Extract_interval(X(j,:)',spike,(delay+t_refr)/dt);
end
b = regress(diff_v_int,XX');                    % linear regression
C = 1/b(3);
g_l = -b(1)*C;
E_l = (b(2)*C)/g_l;
eta = b(4:end)*C;
eta = temp_x'*eta;
time_eta = 0:dt:(length(eta)-1)*dt;
keyboard
%   2. Estimate optimal eta(t) and passive parameters with linear regression
%   3. Estimate moving threshold gamma(t) by maximizing loglikelihood




function [output] = Extract_interval(input,spiketimes,dt)
%
%   Remove the spikes from the data
%   input = vector of real value (i.e. current or voltage)
%   spiketimes = index of the spikes
%   dt = period of time that are removed following each spikes

if(isempty(spiketimes))
    output = input;
else
    output = input(1:spiketimes(1));
    
    for i=1:length(spiketimes)-1
        if(spiketimes(i) + dt <= spiketimes(i+1))
            output = [output;input(spiketimes(i)+dt:spiketimes(i+1))];
        end
    end

    output = [output;input(spiketimes(end)+dt:end)];
end

% OKOKOKOKKOKOKOKOKOKOKOKKOKOKOKOKKOKOOKOKOKOKOKOOKOKOKOKOKOKOOKOKOK %
