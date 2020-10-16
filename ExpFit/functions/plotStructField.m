function h = plotStructField(strcArray,field)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
M=getFieldArray(strcArray,field);
imshow(mat2gray(M))

end

