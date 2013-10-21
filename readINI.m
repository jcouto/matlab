function R = readINI(configfile)
% R = readINI(configfile)
% Reads a configuration (INI) file into a structure R;
% R has the fields of the Nodes of the configuration file
% and each of these is field with the properties.
% Inspired in a piece of code from http://rosettacode.org.
% Note: Ignore empty lines or lines that start with %, # or ;. 
%       Node names are within brackets.
%       Variables follow the notation variable = name
% Joao Couto, October 2013
%

if nargin<1,
    R = [];
    disp('Insuficient inputs to readINI.')
    return;
end;

fid = fopen(configfile);
if fid<0, error('cannot open file %s\n',a); end;
current_node = '';
while ~feof(fid)
    line = strtrim(fgetl(fid));
    % Ignore empty lines or lines that start with %, # or ;. 
    if ~(isempty(line) || all(isspace(line)) || strncmp(line,'#',1)|| strncmp(line,'%',1) || strncmp(line,';',1)),
        [var,tok] = strtok(line,'=');
        var = strtrim(var);
        if(any(var=='['))
            % Node names are within brackets.
            var(find(var == '[' | var == ']')) = [];
            current_node = var;
        end
        % Variables have = simbols.
        if any(tok=='='),
            k = 1;
            while (1)
                [val, tok] = strtok(tok,'');
                val(val == '=' ) = [];
                val = strtrim(val);
                try
                    val = eval(val);
                catch % Handle boolean variables.
                    if strcmp(lower(val),'true');
                        val = logical(1);
                    elseif strcmp(lower(val),'false');
                        val = logical(0);
                    end
                end
                % Store variable in the structure
                if isempty(current_node)
                    R.(var) = val;
                else
                    R.(current_node).(var) = val;
                end
                if isempty(tok), break; end;
                k=k+1;
            end
        end
    end
end

fclose(fid);
