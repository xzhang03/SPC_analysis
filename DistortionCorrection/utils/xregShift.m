function [movingReg, tform] = xregShift(fixed, moving, rtype)
% Aligns the moving image to the fixed image performing only a translation.
% tform is the transformation matrix to align the image.
if nargin < 3
    rtype = 'translation';
end

tform = imregcorr(moving, fixed, rtype);
Rfixed = imref2d(size(fixed));
movingReg = imwarp(moving, tform, 'OutputView', Rfixed);

end
