function [files,kfiles] = list_h5_files(folder, level, sorted)
% [files,kfiles] = list_h5_files(folder, level, sorted)
% Finds all HDF5 files except those that have a *_kernel.dat companion. 
% Returns the files and kernel files sorted by name and in the shape of 
% structures with fields:
%       - path
%       - basename
%       - date (from filename)
%
%       WARNING: Filenames that do not respect the date format yyyymmddHHMMSS
%   will cause datenum to crash.
%
% Example: 
% Returnes all files but the kernels:
%   files = list_h5_files;
% Joao Couto, May 2013

except_folder = {'trash'};
files = {};
kfiles = {};

if ~exist('folder','var');folder = '.';end
if ~exist('level','var');level = 'all';end
if ~exist('sorted','var');sorted = 1;end

% Get the lists of the files
filenames = list_files(folder,'*.h5',level,except_folder,{'*_kernel.dat'});
kfilenames = list_files(folder,'*_kernel*',level,except_folder);
% Get the basename and date from filename.
basenames = cell(size(filenames));
kbasenames = cell(size(kfilenames));
for ii = 1:length(filenames)
    [~, basenames{ii}] = system(sprintf('basename %s',filenames{ii}));
end
for ii = 1:length(kfilenames)
    [~, kbasenames{ii}] = system(sprintf('basename %s',kfilenames{ii}));
end
% Compute the time from the filename
dates = cellfun(@(x)datenum(x,'yyyymmddHHMMSS'),basenames,'UniformOutput',0);
kdates = cellfun(@(x)datenum(x,'yyyymmddHHMMSS'),kbasenames,'UniformOutput',0);
if sorted
    [~,idx] = sort(cell2mat(dates));
    dates = dates(idx);
    filenames = filenames(idx);
    basenames = basenames(idx);
    [~,idx] = sort(cell2mat(kdates));
    kdates = kdates(idx);
    kfilenames = kfilenames(idx);
    kbasenames = kbasenames(idx);
end
if ~isempty(filenames)
    files = cell2struct([filenames,basenames,dates], ...
                        {'path','basename','date'},2);
end
if ~isempty(kfilenames)
    kfiles = cell2struct([kfilenames,kbasenames,kdates], ...
                        {'path','basename','date'},2);
end
end % function