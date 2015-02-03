function printFigWithCaption(figureFile,caption,overlap,rotationAngle)
% printFigWithCaption(figureFile,caption,overwrite file)
% Saves the figure to a file with caption given by string
hideOutput=1; % zero for debug...
% does the file exist

if ~exist(figureFile,'file')
    disp('--> Error [printFigWithCaption]: Can not find figure.')
    disp(['Filename : ',figureFile])
    return
end
[foldername,fname,ext] = fileparts(fullfile(figureFile));
if ~exist('overlap','var')
    overlap = [];
end
if isempty(overlap)
    overlap = false;
end
if ~exist('rotationAngle','var')
    rotationAngle = 0;
end
% get the file absolute path if using ~ to get the user HOME
tmp = strfind(figureFile,'~');
if length(tmp)==1
    if ispc
        figureFile = [getenv('USERPROFILE'),figureFile(tmp+1:end)];
    else
        figureFile = [getenv('HOME'),figureFile(tmp+1:end)];
    end
end
oldPWD = pwd;
cd(foldername)

tmpFname    = 'tempTEX';
% create tmp Latex file.
fid         = fopen([tmpFname,'.tex'],'w');
if ~isempty(caption)
    width = 'width=\textwidth';
else
    width = '';
end
sHeader     = {'\documentclass{article}', ...
    '\usepackage{amsmath}', ...
    '\usepackage[active,tightpage,textmath,displaymath,floats,graphics,previewborder=0.05cm]{preview}', ...
    '\usepackage{graphicx}', ...
    '\begin{document}', ...
    '\begin{figure}', ...
    '\centering', ...
    ['\includegraphics[',width,',angle=',num2str(rotationAngle),']{',fname,ext,'}']};
if ~isempty(caption)
    sHeader = [sHeader,['\caption{',caption,'}']];
end

sHeader = [sHeader,'\end{figure}', ...
    '\end{document}'];
for ii = 1:length(sHeader)
    fprintf(fid,'%s\n',sHeader{ii});
end
fclose(fid);
[tmp,latexPath] = system('which pdflatex');
if isempty(latexPath)
    path1 = getenv('PATH');
    path1 = [path1,':/usr/texbin/']; % On my macOS with MacTex
    setenv('PATH', path1)
end
%pdflatex -jobname here tempTEX.tex
sSubmit = ['pdflatex -interaction=nonstopmode -jobname ',tmpFname,' ',tmpFname,'.tex'];
if hideOutput
    [tmpError,~]=system(sSubmit);
else
    system(sSubmit);
end
%
if exist([tmpFname,'.pdf'],'file')
    if ~overlap
        movefile([tmpFname,'.pdf'],[figureFile(1:end-4),'Caption.pdf'])
    else
        movefile([tmpFname,'.pdf'],figureFile)
    end
    if hideOutput
        delete([tmpFname(1:end-3),'*']) % delete tmp file and other crap that might of been generated...
    end
    
else
    
    disp('printFigWithCaption: Failed usign latex.')
end
cd(oldPWD)

end
