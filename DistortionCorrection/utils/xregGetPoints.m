function [x, y, yFLIM, ySBX] = xregGetPoints(sbx, flimReg, nbpoints)
% Tool to get user defined calibration points to fit a distortion
% correction. The user clicks on 10 landmarks, alternating between the FLIM
% image(left) and the SBX image (right). First click needs to be on the
% FLIM image.

if nargin < 3
    nbpoints = 20;
end

figure
subplot(121)
imshow(flimReg)
title('FLIM')
subplot(122)
imshow(sbx)
title('SBX')
suptitle('Click on corresponding points, starting with the left image')
[points, ~] = ginput(nbpoints);

yFLIM = points(1:2:nbpoints-1);
ySBX = points(2:2:nbpoints);
diffs = yFLIM - ySBX;

y = diffs(:);
x = ySBX(:);
end
