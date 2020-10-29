function [cropped] = crop2reference(im2crop, reference)
% Crops im2crop to be the same size as reference, keeps the center of the
% image
cropped = center_patch(im2crop, size(reference)-[1, 1]);
end
