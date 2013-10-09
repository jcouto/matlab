function [spk,spk_w,t_spk_w,spk_idx,threshold] = extract_extracellular_spikes(data, threshold, time, tpre, tpost, tdead, detection_mode)
% EXTRACT_EXTRACELLULAR_SPIKES Extracts the spikes from extracellular filtered data using the algorithm proposed by Rodrigo Quiroga.
%
% [spk,spkwave,tspkwave,spk_idx,threshold] = extract_extracellular_spikes( DATA, THRESHOLD, T, TPRE, TPOST, TDEAD,detection_mode)
% Channels should be in collums.
%
%  Units:
%
% * DATA               - mV                            []
% * THRESHOLD          - mV                            [-10 mV]
% * TIME VECTOR|DT  - s | s                       [1./30e3 s]
% * TPRE               - ms (time window before spk)   [5ms]
% * TPOST              - ms (time window before spk)   [5ms]
% * TDEAD              - ms (min time between spk)     [2ms]
%
% Note: if threshold is empty, the threshold is selected from the media of
% all points that are larger than 100 mV/ms.
%  Extracts:
% * SPK (Timestamps of the action potentials)
% * SPK_W (Waveforms of the action potentials)%
% * T_SPK_W (Timevector for the waveforms of the action potentials)
% * SPKI_DX (Index of the action potentials)
% Note: SPK and SPK_W are cells if more than one row is provided to
%   DATA.
% Some of this code has been adapted from Quiroga's Waveclus

if ~exist('threshold','var'),thr = [];else thr = threshold;end;
if iscolumn(data)
    data = data';
end
if ~exist('time','var')
    dt = 1./20e3;
    time  = (0:size(data,2)).*dt;
else
    if length(time)>1
        dt = time(2)-time(1);
    else
        dt = time;
        time  = (0:size(data,2)).*dt;
    end
end
if ~exist('tpre','var'),tpre = 1;end
if ~exist('tpost','var'),tpost = 1.5;end
if ~exist('tdead','var'),tdead = 1.5;end
if ~exist('detection_mode','var'),detection_mode = 'neg';end

interpolation = 'y';
int_factor = 2;



N=size(data,1);

spk     = cell(N,1);
spk_w   = cell(N,1);

spk_idx     = cell(N,1);
threshold     = cell(N,1);

ref = ceil(tdead./1000/dt);
w_pre = ceil(tpre./1000/dt);
w_post = ceil(tpost./1000/dt);

% LOCATE SPIKE TIMES
for ii = 1:N
    if isempty(thr)
        thr = find_spike_threshold(data(ii,:));
    end
    threshold{ii} = thr;
    
    switch detection_mode
        case 'pos'
            dd = data;
        case 'neg'
            dd = -data;
        case 'both'
            dd = abs(data);
            
    end
    nspk = 0;
    
    xaux = find(dd(ii,w_pre+2:end-w_post-2) > thr) + w_pre + 1;
    
    xaux0 = 0;
    for i=1:length(xaux)
        if xaux(i) >= xaux0 + ref
            [maxi iaux]=max(dd(ii,xaux(i):xaux(i) + floor(ref/2) - 1));
            nspk = nspk + 1;
            index(nspk) = iaux + xaux(i) -1;
            xaux0 = index(nspk);
        end
    end
    % SPIKE STORING (with or without interpolation)
    ls=w_pre+w_post;
    spikes=zeros(nspk,ls+4);
    xf=data(ii,:);
    xf = [xf zeros(1,w_post)];
    parfor i=1:nspk                          %Eliminates artifacts
        if max(abs( xf(index(i)-w_pre:index(i)+w_post) )) < thr * 50
            spikes(i,:)=xf(index(i)-w_pre-1:index(i)+w_post+2);
        end
    end
    aux = find(spikes(:,w_pre)==0);       %erases indexes that were artifacts
    spikes(aux,:)=[];
    index(aux)=[];
    
    switch interpolation
        case 'n'
            spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation
            spikes(:,1:2)=[];
        case 'y'
            %Does interpolation
            spikes = int_spikes(spikes,w_pre,w_post,detection_mode,int_factor);
    end
    
    spk{ii}     = time(index);
    spk_w{ii}   = spikes;
    
    spk_idx{ii} = index;
end
t_spk_w = linspace(-tpre-dt,tpost,length(spikes(1,:)));
if (N < 2) % then return an array instead of a cell array.
    spk = spk{1};
    spk_w = spk_w{1};
    spk_idx = spk_idx{1};
end


function [spikes1] = int_spikes(spikes,w_pre,w_post,detection_mode,int_factor);
% Interpolates with cubic splines to improve alignment.
% From Quiroga, Rodrigo


ls = w_pre + w_post;

nspk=size(spikes,1);

s=1:size(spikes,2);
ints=1/int_factor:1/int_factor:size(spikes,2);

intspikes=zeros(1,length(ints));
spikes1=zeros(nspk,ls);
switch detection_mode
    case 'pos'
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
            [maxi iaux]=max((intspikes(w_pre*int_factor:w_pre*int_factor+8)));
            iaux = iaux + w_pre*int_factor -1;
            spikes1(i,w_pre:-1:1) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor);
            spikes1(i,w_pre+1:ls) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor);
            
        end
    case 'neg'
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
            [maxi iaux]=min((intspikes(w_pre*int_factor:w_pre*int_factor+8)));
            iaux = iaux + w_pre*int_factor -1;
            spikes1(i,w_pre:-1:1) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor);
            spikes1(i,w_pre+1:ls) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor);
        end
    case 'both'
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
            [maxi iaux]=max(abs(intspikes(w_pre*int_factor:w_pre*int_factor+8)));
            iaux = iaux + w_pre*int_factor -1;
            spikes1(i,w_pre:-1:1) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor);
            spikes1(i,w_pre+1:ls) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor);
        end
end
