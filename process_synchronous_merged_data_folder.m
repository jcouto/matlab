function process_synchronous_merged_data_folder(files, mode)
% process_synchronous_merged_data_folder(files, mode)
% Extracts spike times and voltage traces measured simultaneously.
%    FILES spescifie which files to use (cell array) optional.
%    MODE option is the type of input file:
%        - 'merged_data' :h5 files merged with mcd files using
% the "merge_mcd_to_h5"
% This function does not output any value. Instead it saves a file
% with the same name as the original file with the termination:
%    "_processed_spikes.mat"
% August 2013 - Joao Couto

if ~exist('files','var'),files = {};end
if ~exist('mode','var'),mode = 'merged_data';end

if isempty(files)
    files = list_files('.',sprintf('*.%s',mode),[],{'trash'});
end
file_appendix = '_processed_spikes.mat';
thresh = -20;
tpost = 1.5;
tpre = 1;
for k = 1:length(files)
    fname = files{k};
    [pathstr,filename,~] = fileparts(files{k});
    filename = sprintf('%s%s',filename,file_appendix);
    if ~exist(filename,'file')
        i_process_raw_files(fname,file_appendix,mode)
    end
    load(filename)
    window = int32(0.005/(t(2)-t(1)));
    idx = argfindpeaks(V,thresh,window);
    intra_spk = t(idx);
    intra_pert = findTransitions(I);
    intra_pert = intra_pert(diff(t(intra_pert))<10);
    intra_pert = t(intra_pert(1:2:end));
    if ~exist('raw_psth','var')
        raw_psth = cell(size(spk));
    end
    trig_V = [];
     for ii = 1:length(intra_pert)
         trig_V = [trig_V;V(find(t<=intra_pert(ii)-tpre,1,'last'): ...
             find(t <= intra_pert(ii)+tpost,1,'last'))];
         for jj = 1:length(spk)
            raw_psth{jj} = [raw_psth{jj}, spk{jj}(spk{jj} >= (intra_pert(ii) - tpre) &...
                spk{jj} <= (tpost + intra_pert(ii))) - intra_pert(ii)];
         end
     end
end
save('processed_psth_data.mat','raw_psth','trig_V','intra_pert')
edges = [-tpre*1e3:5:tpost*1e3];
bins = histc([raw_psth{:}]*1e3,edges);
bar(edges,bins,'k')


function i_process_raw_files(fname,file_appendix,mode)

spike_detection_mode = 'both'; %neg, pos, both
tpre = 0.7; % in ms
tpost = 1;
tdead = 1.5;
example_dur = 20; %seconds

if ~exist(fname,'file')
    fprintf(1,'File %s does not exist.\n',fname)
else
    fprintf(1,'Processing file %s',fname);
    switch mode
        case 'merged_data'
            % Open file to read info on the mcd_srate
            [ent,info] = open_h5_trace(fname);
            mcd_srate = info.mcd_srate;
            [ent,info] = load_h5_trace(fname);
            % NOTE: This script is taking the last kernel to make
            %all calculations!
            kfiles = list_files('.','*_kernel.dat');
            Ke = load(kfiles{end});
            V = ent(strcmp({ent.name},'AnalogInput')).data;
            I = zeros(size(V));
            if length(ent) > 1
                I = ent(strcmp({ent.name},'Waveform')).data;
            end
            V = AECoffline(V,I,Ke);
            %V = V - mean(V(I==0));
            t = linspace(0,info.tend,length(V));
            tmcd = linspace(0,info.tend,length(ent(3).data));
            % load data from mcd file
            x = ent(3).data';
            clear('ent')
    end
    xf = zeros(size(x));
    for e = 1:size(x,2)
        xf(:,e) = filter_data(x(:,e),300,5000,mcd_srate);
        % spike detection
        [spk{e},spk_w{e}, t_spk_w,spk_idx{e},...
            threshold{e}] = extract_extracellular_spikes(xf(:,e),...
            [], ...
            tmcd,tpre,...
            tpost,...
            tdead,...
            spike_detection_mode);
        fprintf(1,' . ');
    end
    
    tmp = (tmcd(end) - example_dur)*rand(1);
    idx = (tmcd>tmp &tmcd<tmp+example_dur);
    t_example = tmcd(idx);
    xf_example = xf(idx,:);
    x_example = x(idx,:);
    idx_example = find(idx,1,'first');
    % Data export
    [pathstr,filename,~] = fileparts(fname);
    if isempty(pathstr)
        pathstr = pwd;
    end
    
    save(sprintf('%s/%s%s',pathstr,filename,file_appendix),'t','V',...
        'I','spk','t_spk_w','spk_w','spk_idx','threshold','x_example',...
        'xf_example','t_example','idx_example')
    clear('t','V','x','xf',...
        'I','spk','t_spk_w','spk_w','spk_idx','threshold','x_example',...
        'xf_example','t_example','idx_example')
    fprintf(1,' Saved %s%s file.\n',filename,file_appendix);
end