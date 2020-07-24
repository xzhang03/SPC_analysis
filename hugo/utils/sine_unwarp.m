function [out_image,D] = sine_unwarp(image,a,w)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
maxI = max(image(:));

yvect = 1:size(image,2);
proj = a * sin(1/w * (yvect - size(image,2)/2));

D = zeros(size(image,1),size(image,2),2);
D(:,:,1) = repmat(proj,[size(image,1),1]);

out_image = imwarp(image,D,'cubic');
out_image(out_image>maxI) = maxI;
end

