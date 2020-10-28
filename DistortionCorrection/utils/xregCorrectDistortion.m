function [flimCorrected] = xregCorrectDistortion(flimReg, D)
% Applies a pre-computed displacement field to an image to correct for
% scanning induced distortions.

flimCorrected = imwarp(flimReg, D, 'cubic');

end
