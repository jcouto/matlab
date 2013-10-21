function [outvar] = process_eif_folder(foldername)

outvar = [];

if ~exist('foldername','var')
    foldername = cd(cd('./'));
end
    files = list_h5_files(foldername);

for ii = 1:length(files)
    process_eif_file(files(ii).path)
end
