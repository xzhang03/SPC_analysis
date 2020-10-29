function [tmCorrected] = xregCorrectTM(tm, tform, D)
% Applies the shift (tform) and distortion correction (D) to a tm image.
Rfixed = imref2d(size(tm));
tmShifted = imwarp(tm, tform, 'OutputView', Rfixed);
tmCorrected = xregCorrectDistortion(tmShifted, D);

end
