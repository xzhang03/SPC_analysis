function [im1cropped, im2cropped] = xregCropping(im1, im2, margin)
% Crops two images to the same size using the margin defined by the user to
% get rid of deadbands on edges of the image.
if nargin < 3
    margin = 100;
end

w = [min(size(im1, 1), size(im2, 1)) - margin, min(size(im1, 2), ...
    size(im2, 2)) - margin];
im1cropped = center_patch(im1, w);
im2cropped = center_patch(im2, w);
end
