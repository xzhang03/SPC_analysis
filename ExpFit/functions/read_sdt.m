function [photCount,outImage] = read_sdt(path,factor,mode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if mode == 'matfile'
    [fileDir,fileName,fileExt] = fileparts(path);
    matpath = fullfile(fileDir,strcat(fileName,'.mat'));
    outImage = load(matpath);
    outImage = outImage.data;
elseif mode == 'sdtfile'
    sdt = bfopen(path);
    outImage = cat(3,sdt{1}{:,1});
end

if factor
    outImage = bin2d(outImage,factor);
end

photCount = squeeze(sum(sum(outImage,1),2));

end

