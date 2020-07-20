function [np, np_size] = MakeNPRing(im, inner, outter, im_neg, im_allow, minpixels)
% MakeNPRing makes neurpil ring given a inner distance, a starting outter
% distance, a negative image to subtract from, an image of allowable area, 
% and a minimal number of pixels.
% [np, np_size] = MakeNPRing(im, inner, outter, im_neg, im_allow, minpixels)

% Initialize
np_size = 0;
outter = outter - 1;

% Binarize
im = im > 0;
im_neg = im_neg > 0;

% Inner strel
in = strel('disk', inner);
im_in = imdilate(im, in);

while np_size < minpixels
    % Propagate outter distance
    outter = outter + 1;
    out = strel('disk', outter);
    
    % Candidate neuropil
    np = ((imdilate(im, out) - im_in - im_neg) .* im_allow) > 0;
    np_size = sum(np(:));
    
end

end