function [x, y, yFLIM, ySBX] = xregGetPoints(sbx, flimReg, nbpoints)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
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
