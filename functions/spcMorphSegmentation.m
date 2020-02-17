function spcMorphSegmentation(mouse, date, run, varargin)
% spcMorphSegmentation uses morphological filtering to segment ROIs

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Use registered?
addOptional(p, 'useregistered', true);

% Binning
addOptional(p, 'binxy', 1); % Binning

% Segmentation
addOptional(p, 'Threshold', 0.05); % Threshold
addOptional(p, 'SizeThresh', [30 300]); % Size threshold (after binning)

% Take a look
addOptional(p, 'viewresults', true);


% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Case and type
mouse = upper(mouse);
if ~ischar(date)
    date = num2str(date);
end

% User (add yourself if needed)
if isempty(p.user)
    switch mouse(1:2)
        case 'SZ'
            p.user = 'stephen';
        case 'AL'
            p.user = 'andrew';
        case 'HK'
            p.user = 'hakan';
        case 'YL'
            p.user = 'yoav';
            
    end
end

%% IO
% Get paths
spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);

% Load
if p.useregistered
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.regtif_photons));
else
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_photons));
end

%% Segmentation
% Bin if needed
if p.binxy > 1
    im_photon = binxy(im_photon, p.binxy);
end

% Frame to segment
im2seg = median(im_photon, 3);

% Local normalize
f_prime = im2seg - imgaussfilt(im2seg, 8);
im2seg2 = f_prime ./ (imgaussfilt(f_prime.^2, 30).^(1/2));

% Segmentation
icaguidata = sbxMorphologicalFilterExtractROIs_SZ(im2seg2, 'PlotOrNot', false,...
    'Threshold', p.Threshold, 'SizeThresh', p.SizeThresh);

% Take a look
if p.viewresults
    rgb = repmat(mat2gray(im2seg2), [1 1 3]);
    rgb(:,:,3) = 0;
    rgb(:,:,2) = edge(icaguidata.AllFilters > 0);
    figure;imshow(rgb)
end

%% Save
save(fullfile(spcpaths.fp_out, spcpaths.mat), 'im2seg', 'im2seg2', 'icaguidata');
end
