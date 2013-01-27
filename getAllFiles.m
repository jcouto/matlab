function [fileList,folderList] = getAllFiles(dirName)
    % Usage:
    % [fileList,folderList] = getAllFiles(dirName | [pwd])
    % Lists all files and folders under a specified path.
    % With no input: returns all files under current working directory
    % Filters are NOT IMPLEMENTED yet...
    
  if ~exist('dirName','var')
      dirName = pwd;
  end
  dirData = dir(dirName);                           % Get the data for the current directory
  dirIndex = [dirData.isdir];                       % Find the index for directories
  fileList = {dirData(~dirIndex).name}';            % Get a list of the files
  
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),... % Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};               % Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});       % Find index of subdirectories
                                                    %   that are not '.' or '..'
  folderList = cellfun(@(x)fullfile(dirName,x),...
      subDirs(find(validIndex)),'uniformoutput',false)';
  
      
  for iDir = find(validIndex)                       % Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});      % Get the subdirectory path
    [tmpfiles, tmpfolders] = getAllFiles(nextDir);
    fileList = [fileList; tmpfiles];                % Recursively call getAllFiles
    folderList = [folderList; tmpfolders];
  end
end