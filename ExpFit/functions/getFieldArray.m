function M = getFieldArray(structArray,field)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
structArray=squeeze(structArray);
sizeImage = size(structArray);
M = reshape([structArray.(field)],sizeImage(1),sizeImage(2));
end

