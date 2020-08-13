function [movingReg,tform] = xregShift(fixed,moving,rtype)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    rtype = 'translation';
end

tform = imregcorr(moving,fixed,rtype);
Rfixed = imref2d(size(fixed));
movingReg = imwarp(moving,tform,'OutputView',Rfixed);

end

