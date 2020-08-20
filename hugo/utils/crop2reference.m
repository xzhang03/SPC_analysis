function [cropped] = crop2reference(im2crop,reference)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cropped = center_patch(im2crop,size(reference)-[1,1]);
end

