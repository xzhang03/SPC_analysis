function [photCount,outImage] = read_sdt(path,factor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    factor = false;
end

sdt = bfopen(path);
outImage = cat(3,sdt{1}{:,1});

if factor
    outImage = bin2d(outImage,factor);
end

photCount = squeeze(sum(sum(outImage,1),2));

end

