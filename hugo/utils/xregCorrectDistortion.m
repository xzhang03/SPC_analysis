function [flimCorrected] = xregCorrectDistortion(flimReg, D)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

flimCorrected = imwarp(flimReg,D,'cubic');

end

