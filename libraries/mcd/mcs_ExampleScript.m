% mcs_ExampleScript  
% nsMCDlibrary.dll Neuroshare implementation for MC_Rack data 
%
% mcs_plot generates a graphical output and gives additional information about the data 
% mcs_plot produces all parameters to call and perfom mcs_GetEntities ,mcs_Graphic and mcs_Info
% 
% Usage:
%       /
%
% Description: mcs_plot produces parameters to call and perform
%              mcs_GetEntities, mcs_grafik and mcs_Info
%
%
% Information you have to give: name of DLL "*.dll"
%                               name or pathway of MC_Rack data file "*.mcd"                       
%
% Parameters produced:
%               hFile	        Handle/Indentification number to an open file.
%               FileInfo            
%               LibraryInfo
%
% Return Values:
%               /
                        

close all
clear all
clc

itemCount = 250000;                                                        % set number of itemCount (e.g.250000)
%-------------------------------------------------------------------------------------------------------------------------------------------------------

[ns_RESULT] = mcs_SetLibrary('nsMCDlibrary.dll');                          % first step: Opens a Neuroshare shared Library, load the appropriate DLL

if (ns_RESULT ~= 0)                                                        % error message if nsresult is unequal 0
    disp('DLL was not found!');
    return
end

[ns_RESULT, LibraryInfo] = ns_GetLibraryInfo();                            % Obtains information about the library
 

if (ns_RESULT ~= 0)
    disp('Library information did not load!');                             % error message if nsresult is unequal 0
    return
end



[ns_RESULT, hfile] = ns_OpenFile('z:\MC_Rack Data\mailin.mcd');

%[ns_RESULT, hfile] = ns_OpenFile('Z:\MC_Rack Data\DataCardio.mcd');       % second step: Opens a neural data file, specified by filename
                                                                           %(e.g.'Z:\MC_Rack Data\DataCardio.mcd' via pathway) 
                                                                           % and returns a file handle called hfile, that contains an identification number 
if (ns_RESULT ~= 0)
    disp('Data file did not open!');                                       % error message if nsresult is unequal 0
    return
end


[ns_RESULT,FileInfo]=ns_GetFileInfo(hfile);                                % ns_GetFileInfo uses given "hfile" to generate "FileInfo", that contains some information about the file 
                                                                           % like FileType, EntityCount, TimeStampResolution, TimeSpan, AppName,Time_Year, 
                                                                           % Time_Month, Time_Day, Time_Hour, Time_Min, Time_Sec, Time_MilliSec, FileComment 
if (ns_RESULT ~= 0)
    disp('Data file information did not load!');                           % error message if nsresult is unequal 0
    return
end

% --------------------------------------------------------------------------------------------------------------------------------------------------                                                                   
LibraryInfo

FileInfo

hfile                                                                      % output is given in CommandWindow 

 %--------------------------------------------------------------------------------------------------------------------------------------------------

entities = mcs_GetEntities(hfile, 'elec0001');                             % access to file "mcs_GetEntities", using "hfile" and 
                                                                           % stream name (e.g.'elec0001'addressses electrode raw data) returns "entities"

% mcs_grafik(hfile,entities,0,100000)                                      % access to file "mcs_grafik", using "hfile", "entities", firstItem (e.g.0), itemCount(e.g.100000) creates graphic 
[ns_RESULT,entityInfo]= ns_GetEntityInfo(hfile,entities);

entityInfo(1).ItemCount                                                    % ItemCount needed to call mcs_Graphic

mcs_Graphic(hfile,entities,1,entityInfo(1).ItemCount)                      % access to file "mcs_Graphic", using "hfile", "entities", firstItem (e.g.0), itemCount(e.g.100000) creates graphic

b = mcs_Info(hfile);                                                       % access to file "mcs_Info", using "hfile", first entry of varargout is defined as b 

%---------------------------------------------------------------------------------------------------------------------------------------------------------

b                                                                          % b is shown via CommandWindow, containing matrix m

%----------------------------------------------------------------------------------------------------------------------------------------------------------

