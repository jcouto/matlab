function [data] = load_trace_from_mcd_file(hfile,ChannelType,channel,startend,loadByIndex)
%[data] = load_trace_from_mcd_file(hfile,ChannelType,channel,startend)
% Loads the data from an mcd file.
if ~exist('ChannelType','var')
    ChannelType = 'elec';
    %ChannelType = 'digi';
end

% Get file information

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
%{EntityInfo(AnalogList).EntityLabel}'
tmp = cellfun(@(x)x(1:4),{EntityInfo(AnalogList).EntityLabel},'uniformoutput',0);
ChannelList = AnalogList(find(strcmp(tmp,ChannelType)));
ChannelNumber = cellfun(@(x)str2num(x(14:18)),{EntityInfo(ChannelList).EntityLabel})';
channel_to_load = find(ChannelNumber==channel);
if isempty(channel_to_load)
    sprintf('Could not find channel %d. Channels available are: ',channel)
    data = ChannelList;
    return
end
% How many of a particular entity do we have
% cNeural = length(NeuralList)
% cSegment = length(SegmentList)
% cAnalog = length(AnalogList)
% cEvent = length(EventList)

% Load the data from the selected channel
EntityID = ChannelList(channel_to_load);
[nsresult, analogInfo] = ns_GetAnalogInfo(hfile, EntityID);

istartend = nan(2,1);

if ~exist('loadByIndex','var')
    [ns_result, istartend(1)] = ns_GetIndexByTime(hfile, EntityID, startend(1),0);
    [ns_result, istartend(2)] = ns_GetIndexByTime(hfile, EntityID, startend(2),0);
else
    istartend = startend;
end

[nsresult,cont_count, x] = ns_GetAnalogData(hfile,EntityID,istartend(1),diff(istartend));

data.entityID = EntityID;
data.data = (x./(analogInfo.Resolution*1e3));
data.xpos = analogInfo.LocationX;
data.ypos = analogInfo.LocationY;
data.zpos = analogInfo.LocationZ;
data.startend = startend;
data.dt = 1./analogInfo.SampleRate;
data.channel = channel;
data.date = date;
