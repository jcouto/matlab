function [hfile,outinfo] = open_mcd_file(filename)
% Opens a multichannel systems file. Use load_trace_from_mcd_file(file,channel,startend)
% to load a channel.
if ~ischar(filename)
   errormsg = sprintf('Argument "filename" needs to be a string.');
   error(errormsg)
end
if ~exist(filename, 'file') == 2
   errormsg = sprintf('File %s not found!',filename);
   error(errormsg)
end

open_neuroshare()
% Open file
[nsresult, hfile] = ns_OpenFile(filename);
if (nsresult ~= 0)
    disp('Data file did not open!');
    hfile = [];
    return
end
% Return information
[nsresult, info] = ns_GetFileInfo(hfile);
% Gives you EntityCount, TimeStampResolution and TimeSpan
if (nsresult ~= 0)
    disp('Data file information did not load!');
    data = [];
    return
end

%
date = datestr(datenum(info.Time_Year,info.Time_Month,...
    info.Time_Day,info.Time_Hour,...
    info.Time_Min,info.Time_Sec+info.Time_MilliSec/1e3));
tend = info.TimeSpan; % seconds
if ~exist('startend','var')
    startend = [0,tend];
end

% Build catalogue of entities
[nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : info.EntityCount]);

NeuralList = find([EntityInfo.EntityType] == 4);    % List of EntityIDs needed to retrieve the information and data
SegmentList = find([EntityInfo.EntityType] == 3);
AnalogList = find([EntityInfo.EntityType] == 2);
EventList = find([EntityInfo.EntityType] == 1);
tmp = cellfun(@(x)x(1:4),{EntityInfo(AnalogList).EntityLabel},'uniformoutput',0);
% Split electrode raw data and digital channels
ElectrodeList = AnalogList(find(strcmp(tmp,'elec')));
DigitalList = AnalogList(find(strcmp(tmp,'digi')));
% Retrieve channel indexes.
outinfo.electrodes = cellfun(@(x)str2num(x(15:18)),{EntityInfo(ElectrodeList).EntityLabel},'uniformoutput',1);
outinfo.digital = cellfun(@(x)str2num(x(15:18)),{EntityInfo(DigitalList).EntityLabel},'uniformoutput',1);
[ns_result, outinfo.nsamples] = ns_GetIndexByTime(hfile, ElectrodeList(1), tend+100,0);

outinfo.date = date;
outinfo.filename = filename;
outinfo.tend = tend;
outinfo.dt = info.TimeStampResolution;


function open_neuroshare()
% Imports the neuroshare lib.
    switch computer()
        case 'PCWIN'
            [pathname, name, ext]=fileparts(which('nsMCDLibrary.dll'));
         case 'PCWIN64'
            [pathname, name, ext]=fileparts(which('nsMCDLibrary64.dll'));
        case 'GLNX86'
            [pathname, name, ext]=fileparts(which('nsMCDLibraryLinux32.so'));
        case 'GLNXA64'
            % TODO: Will need appropriate mexprog
            [pathname, name, ext]=fileparts(which('nsMCDLibraryLinux64.so'));
        case {'MACI', 'MACI64'}
            [pathname, name, ext]=fileparts(which('nsMCDLibrary.dylib'));
    end
    if (ns_SetLibrary([pathname filesep name ext]) ~= 0)
        error('''%s'' was not found on the MATLAB path',...
            fullfile(pathname, name, ext));
        return
    end
