function h = plotStructField(strcArray,field)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
strcArray=squeeze(strcArray);
sizeImage = size(strcArray);
M = reshape([strcArray.(field)],sizeImage(1),sizeImage(2));
imshow(mat2gray(M))

end

