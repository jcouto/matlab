function [var] = process_eif_folder(foldername,C,FORCE_FROM_SCRATCH,threshold,plotvar)
% Processes [var] = process_eif_folder(foldername,C,FORCE_FROM_SCRATCH,plotvar)

if ~exist('C','var');C = []; end
if ~exist('FORCE_FROM_SCRATCH','var');FORCE_FROM_SCRATCH = 0; end
if ~exist('plotvar','var');plotvar = 1; end
if ~exist('threshold','var');threshold = []; end

if ~exist('foldername','var')
    foldername = cd(cd('./'));
end
files = list_h5_files(foldername);

for ii = 1:length(files)
    [folder,fname] = fileparts(files(ii).path);
    tmp = list_files(folder,sprintf('%s*eif.mat',fname));
    if isempty(tmp) | FORCE_FROM_SCRATCH
        tmp = process_eif_file(files(ii).path,[],C,threshold,plotvar);
    else
        fprintf(1,'Skipping file %s\n',files(ii).path)
    end
end
%% Get experiment name
expName = regexp(pwd,'[0-9]{8}[A-Z][0-9]{2}','match');
if isempty(expName)
    try
        % Try to get the name from two folders above.
        [~,expName] = fileparts(fileparts(fileparts(pwd)));
    catch
        expName = {'unknown'};
    end
end
if iscell(expName)
   expName = expName{end};
end
%tmp = strsplit(folder,expName);
tmp = regexp(folder,expName,'split');

tmp = tmp{end};
tmp(tmp =='/') = '_';
appendix = sprintf('%s%s',expName,tmp);


%% Merge all files
files = list_files(foldername,'*_eif.mat','all',{'trash'});
if ~isempty(files)
    for ii = 1:length(files)
        data(ii) = load(files{ii});
    end
    
    tmp = [data.capacitance];
    var.C = [tmp.C];
    var.Vrest = [tmp.Vrest];
    tmp = [data.eif];
    var.param = vertcat(tmp.param);
    tmp = [data.reif];
    %%
    ref_t = [];
    ref_param = [];
    for ii = 1:length(tmp)
        ref_t = [ref_t, tmp(ii).windows(1:end-1)];
        %     ref_t = [ref_t, tmp(ii).windows(1:end-1)+diff(tmp(ii).windows)/2];
        ref_param = vertcat(ref_param,tmp(ii).param);
    end
    %%
    var.ref_param = ref_param;
    var.ref_t = ref_t;
    
    if plotvar
        
        fig = figure('visible','on');clf
        ax(1) = subplot(2,5,1); % Capacitance
        edges = [0:10:500];
        [count] = histc(var.C,edges);
        tmp = bar(edges(count>0),count(count>0));
        set(tmp,'facecolor','k','facecolor','k')
        axis tight
        tmp = xlim;
        xlim([min(tmp)-diff(tmp)/2,max(tmp)+diff(tmp)])
        xlabel('Capacitance (pF)')
        ylabel('Counts')
        
        label = {'Time constant (\tau_m)','Resting potential (E_m)', 'Spike width (/delta_T)','Spike threshold (V_t)'};
        par_edges = {[0:3:100],[-90:3:-30],[0:0.3:20],[-90:3:-30]};
        
        %%
        for ii = 1:size(var.param,2) % eif parameters
            ax(ii+1) = subplot(2,5,ii+1);
            [count] = histc(var.param(:,ii),par_edges{ii});
            tmp = bar(par_edges{ii}(count>0),count(count>0));
            set(tmp,'facecolor','k','facecolor','k')
            axis tight
            tmp = xlim;
            xlim([min(tmp)-diff(tmp)/2,max(tmp)+diff(tmp)])
            xlabel(label{ii})
            ylabel('Counts')
            
        end
        %%
        ax(2+size(var.param),2) = subplot(2,5,ii+2);
        edges = [-100:3:30];
        [count] = histc(var.Vrest,edges);
        tmp = bar(edges(count>0),count(count>0));
        set(tmp,'facecolor','k','facecolor','k')
        axis tight
        tmp = xlim;
        xlim([min(tmp)-diff(tmp)/2,max(tmp)+diff(tmp)])
        xlabel('V_{rest} (mV)')
        ylabel('Counts')
        label = {'Time constant (\tau_m)','Resting potential (E_m)',...
            'Spike width (/delta_T)','Spike threshold (V_t)'};
        ref_f = {@(tau,p,v)1.0/tau + p(1) * exp(-v/p(2)),...
            @(Em,p,v)Em + p(1) * exp(-v/p(3)) + p(2) * exp(-v/p(4)),...
            @(VT,p,v)VT + p(1) * exp(-v/p(2)),...
            @(dT,p,v)dT + p(1) * exp(-v/p(2))};
        x0 = {[0.1,20],[-30,30,20,50],[mean(var.param(:,3)),1e-5],[10,20]};
        lb = {[0,0.1],[-100,-100,1e-5,1e-5],[-20,1e-3],[0,0.1]};
        ub = {[3,300],[100,100,250,250],[+20,200],[100,300]};
        options = optimset('TolFun',1e-6,'maxiter',100,'display','off');
        for ii = 1:size(var.param,2)
            ax(ii+2+size(var.param,2)) = subplot(2,5,ii+2+size(var.param,2));
            X = var.ref_t;
            Y = var.ref_param(:,ii)';
            if (ii == 1) % invert tau
                Y = 1./Y;
            end
            [X,idx] = sort(X);
            Y = Y(idx);
            plot(X,Y,'ko','markersize',3)
            hold on
            % figure(9),clf
            %  plot(X,-Y/C,'ko')

            [ref_x{ii},resnorm,r,~,o] = lsqcurvefit(@(p,v)ref_f{ii}(mean(var.param(:,ii)),p,v),x0{ii},X,Y,lb{ii},ub{ii},options);
            
            
            
            edges = linspace(0,250,100);
            plot(edges,ref_f{ii}(mean(var.param(:,ii)),ref_x{ii},edges),'r')
            %         [~,edges] = max(arrayfun(@(x)length(x.windows),[data.reif]));
            %         edges = [data(edges).reif.windows,max(data(edges).reif.windows)+50];
            %         [x,y,ys] = binSamples(var.ref_t,var.ref_param(:,ii),edges);
            %         errorbar(x,y,ys)
            
            ylabel(label{ii})
            xlabel('time from spike (ms)')
            
        end
        var.ref_x = ref_x;
        var.ref_f = cellfun(@(x)func2str(x),ref_f,'uniformoutput',0);
        appendix = 'refractory_param_summary.pdf'
        figName = sprintf('%s/%s.pdf',foldername,appendix);
        print(gcf,'-dpdf',figName)
        
    end
    
else
    var = [];
end
%%

save(sprintf('%s_reif_param.mat',expName),'-struct','var');

