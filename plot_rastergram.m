function out=plot_rastergram(sp,offset,increment,type,varargin)
% PLOT_RASTERGRAM Plot a rastergram from spike data.
%
% PLOT_RASTERGRAM(SP) with sp can be a binary array or a cell array, in any case the spike times
% corresponds to a vertical line; the rows correspond to different trials or cells.
%
%
if isempty(offset)
    offset = 0;
end
if isempty(increment)
    increment = 1;
end
if isempty(type)
    type = 'convimg';
end
if ~iscell(sp)
    if ~(length(unique(sp(:))) < 4)
        % Then it is probably not a  a binary spiketrain...
        error('First input must me a cell array or a binary matrix of spiketimes.')
    end
end

switch type
    case 'convimg'
        if iscell(sp)
            sp = binary_spiketrains(sp,[],0.001);
        end
        cidx = find(strcmp(varargin,'color'));
        if ~isempty(cidx)
            cc = 1 - varargin{cidx+1};
        else
            cc = [1,1,1];
        end
        if size(cc,1) == 1
            cc = repmat(cc,size(sp,1),1);
        end
        %%
        x = -100:100;
        sigma = 3;
        kernel = exp(-x .^ 2 / (2*sigma ^ 2)); kernel = kernel./sum(kernel);
%         kernel(x>-2 & x<3) = 1; % tap the top
        for i = 1:size(sp,1)
            sp(i,:) = filtfilt(kernel,1,sp(i,:));
        end
        
        sp(sp>0.2) = 0.2;
        
        sp = sp./max(sp(:));
        im = cat(3,1-sp.*repmat(cc(:,1),1,size(sp,2)),1-sp.*repmat(cc(:,2),1,size(sp,2)),1-sp.*repmat(cc(:,3),1,size(sp,2)));
        try
        imagesc([0,size(sp,2)-1],offset+[0,size(sp,1)-1].*increment,im)
        catch
            warning('A raster failed...')
        end
    case 'blackimg'
        if iscell(sp)
            sp = binary_spiketrains(sp,[],0.001);
        end
        cidx = find(strcmp(varargin,'color'));
        if ~isempty(cidx)
            cc = 1 - varargin{cidx+1};
        else
            cc = [0,0,0];
        end
        x = -100:100;
        sigma = 3;
        kernel = exp(-x .^ 2 / (2*sigma ^ 2)); %kernel = kernel./sum(kernel);
        kernel(x>-2 & x<3) = 1; % tap the top
        for i = 1:size(sp,1)
            sp(i,:) = filtfilt(kernel,1,sp(i,:));
        end
        sp = sp;
        colormap(1-gray)
        imagesc(sp)
    case 'line'
        if ~iscell(sp)
            
            minimum=-0.5;
            maximum=0.5;
            if (size(cc,1) ~= length(spikes))
                cc = repmat(cc,length(spikes),1);
            end
            for ii=1:length(spikes)
                for jj=1:length(spikes{ii})
                    line(spikes{ii}(jj).*[1,1],[minimum,maximum]+counter,'color',cc(ii,:),'linewidth',linesize)
                end
                counter=counter+1;
            end
        end
    otherwise
        disp('---> What is it that you want? Options are (massive) and (few)...')
        return
end
%xlabel('Time (s)','fontname','Arial','fontsize',12,'fontweight','bold')
%ylabel('Trials','fontname','Arial','fontsize',12,'fontweight','bold')