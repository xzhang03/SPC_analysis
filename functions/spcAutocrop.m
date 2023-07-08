function bbox = spcAutocrop(im, varargin)

%% Parse inputs
if nargin < 2
    varargin = {};
end

p = inputParser;
addOptional(p, 'size', 3);
addOptional(p, 'threshold', 1);
addOptional(p, 'buffer', [10 10 10 10]); %[col_start, row_start, col_end, row_end]

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Sizes
% Size
imsize = size(im);

% Threshold
imthresh = imclose(im,strel('disk',p.size)) >= p.threshold;

% Label and area
imlabel = bwlabel(imthresh, 4);
props = regionprops(imlabel, 'Area', 'BoundingBox');
areas = [props(:).Area];

% Keep
[~, maxind] = max(areas);
bbox = round(props(maxind).BoundingBox);

% Buffer
bbox(3:4) = bbox(1:2) + bbox(3:4) + p.buffer(1:2);
bbox(1:2) = bbox(1:2) - p.buffer(1:2);

% Check out of bounds
bbox(1) = max(bbox(1), 1);
bbox(2) = max(bbox(2), 1);
bbox(3) = min(bbox(3), imsize(2));
bbox(4) = min(bbox(4), imsize(1));

end