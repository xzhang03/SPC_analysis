function [tmCorrected] = xregCorrectTM(tm, tform, D)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Rfixed = imref2d(size(tm));
tmShifted = imwarp(tm, tform, 'OutputView', Rfixed);
tmCorrected = xregCorrectDistortion(tmShifted, D);

end
