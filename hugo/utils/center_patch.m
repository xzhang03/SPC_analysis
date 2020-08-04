function [P] = center_patch(I,w)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if numel(w)==1
    xrange = round(size(I,1)/2 - w/2) : round(size(I,1)/2 + w/2);
    yrange = round(size(I,2)/2 - w/2) : round(size(I,2)/2 + w/2);
else 
    xrange = round(size(I,1)/2 - w(1)/2) : round(size(I,1)/2 + w(1)/2);
    yrange = round(size(I,2)/2 - w(2)/2) : round(size(I,2)/2 + w(2)/2);
end
P = I(xrange,yrange);
end

