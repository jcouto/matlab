function R = readINI(configfile)
% R = readINI(configfile)
% Reads a configuration (INI) file into a structure R;
% R has the fields of the Nodes of the configuration file
% and each of these is field with the properties.
% Inspired in READCONF from Matlab Central. 
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
    if ~(isempty(line) || all(isspace(line)) || strncmp(line,'#',1) || strncmp(line,';',1)),
	[var,tok] = strtok(line,'=');
    var = strtrim(var);
    if(any(var=='['))
%         fprintf(1,'Found node %s',var)
        var(find(var == '[' | var == ']')) = [];
        current_node = var;
    end
	if any(tok=='='),
		k = 1; 
		while (1)
			[val, tok] = strtok(tok,'');
            val(val == '=' ) = [];
            val = strtrim(val);
            try
                val = eval(val);
            catch
                if strcmp(lower(val),'true');
                    val = logical(1);
                elseif strcmp(lower(val),'false');
                    val = logical(0);
                end
            end
            if isempty(current_node)
                R.(var) = val;  	% return value of function 
            else
                R.(current_node).(var) = val;  	% return value of function 
            end
		if isempty(tok), break; end;
			k=k+1;
		end;
	end;
    end;
end; 
fclose(fid);
