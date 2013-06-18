function process_sinusoid_modulation_file(file)
% This function processes a file that had a sinusoidal waveform injected.

[files,kfiles] = list_h5_files;

if ~exist('plotvar','var');plotvar = 1;end
if ~exist('ii','var');ii = length(files);end

[tall,Vall,Iall,metadata,info] = i_load_compensated_voltage(files(ii),kfiles);


% Assumes the last thing is the sinusoids protocol (and that there is a
% tail just before the end).
tprot = unique(cumsum(metadata(:,1)));
tstart_sin = tprot(end-2);
tend_sin = tprot(end-1);
V0 = mean(Vall(tall<tprot(1)))

dt = tall(2)-tall(1);
dVdt = (diff(Vall)*1.0e-3)./dt;
[~,~,~, spkidx] = extract_spikes( dVdt, 100, tall, 1, 3, 3);

keyboard
metadata

function [t, V, I, metadata, info] = i_load_compensated_voltage(file,kfiles)
% Internal function to load voltage trace
V = [];
I = [];
metadata = [];
[ent, info] = load_h5_trace(file.path);
idx = find(strcmp('RealNeuron',{ent.name}));
idx = [idx, find(strcmp('AnalogInput',{ent.name}))];
V = [ent(idx(1)).data];
t = linspace(0,info.tend,length(V));
% Do we need AEC?
if isempty(ent(idx(1)).metadata) && ~isempty(kfiles)
    % Is there I?
    idx = find(strcmp('Waveform',{ent.name}));
    metadata = ent(idx).metadata;
    % if there was a holding potential include it on the AEC current
    idx = [idx,find(strcmp('Constant',{ent.name}))];
    I = sum(vertcat(ent(idx).data),1);
    [~,k]  = min((file.date) - [kfiles.date]);
    Ke=load(kfiles(k).path);
    V = AECoffline(V,I,Ke);
end