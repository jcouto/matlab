function importlib(name,helper)
% imports particular libraries, the location of these libraries are
% specified in the function file. - edit importlib
% Available libraries:
%       - 'prct'
%       - 'csprc'

AVAILABLE = {'prct','csprc'};
PATHTOLIB = {{'/Users/joao/projectPurkinjePrc/matlablib'},...
    {'/Users/joao/projectPurkinjePrc/csprc/',...
    '/Users/joao/projectPurkinjePrc/csprc/functions',...
    '/Users/joao/projectPurkinjePrc/csprc/functions/utils',...
    '/Users/joao/projectPurkinjePrc/csprc/functions/utils/spikes',...
    '/Users/joao/projectPurkinjePrc/csprc/functions/model_selection'}};
flag = 0;

if ~exist('name','var')
    print_help(AVAILABLE);
    return
else
    idx = find(cellfun(@(x)strcmp(x,name),AVAILABLE));
    if ~isempty(idx)
        flag = 1;
        for pp = PATHTOLIB{idx}
            path(path,pp{1})
            fprintf(1,'Added %s to path.\n',pp{1})
            if exist('helper','var')
                help(pp{1})
            end
        end
        
        
    end
end
if ~flag
    fprintf(1,'Library %s not found.\n',name)
    print_help(AVAILABLE)
end
end


function print_help(AVAILABLE) 
    fprintf(1, '\nAvailable libraries:\n')
    for ll = AVAILABLE 
        fprintf(1, '\t\t\t - %s\n',ll{1})
    end
end