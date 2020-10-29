function [P] = center_patch(I, w)
% Removes w pixels from all 4 edges of image I. w can also be a vector to
% remove bands of different sizes in x and y directions.
if numel(w) == 1
    xrange = floor(size(I, 1)/2-w/2):floor(size(I, 1)/2+w/2);
    yrange = floor(size(I, 2)/2-w/2):floor(size(I, 2)/2+w/2);
else
    xrange = floor(size(I, 1)/2-w(1)/2):floor(size(I, 1)/2+w(1)/2);
    yrange = floor(size(I, 2)/2-w(2)/2):floor(size(I, 2)/2+w(2)/2);
end
P = I(xrange, yrange, :);
end
