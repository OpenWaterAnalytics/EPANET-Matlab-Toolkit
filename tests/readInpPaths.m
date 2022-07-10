function [InpList] = readInpPaths(DataFolder)
% This function is based on
% https://www.mathworks.com/matlabcentral/fileexchange/31343-enlist-all-file-names-in-a-folder-and-it-s-subfolders
% Author: Thokare Nitin D.
% Modified: Marios S. Kyriakou
%
% Example:
% [InpList] = readInpPaths(pwd)

DirContents=dir(DataFolder);
InpList=[];

for i=1:numel(DirContents)
    if(~(strcmpi(DirContents(i).name,'.') || strcmpi(DirContents(i).name,'..')))
        if(~DirContents(i).isdir)
            extension=DirContents(i).name(end-2:end);
            if(numel(find(strcmpi(extension, 'inp')))~=0)
                InpList=cat(1,InpList,{[DataFolder,'/',DirContents(i).name]});
            end
        else
            getlist=readInpPaths([DataFolder,'/',DirContents(i).name]);
            InpList=cat(1,InpList,getlist);
        end
    end
end
end
