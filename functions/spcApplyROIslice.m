function spcApplyROIslice(mouse, date, run, varargin)
% spcApplyROI applies ROIs to FLIM images. The slice version generates a
% matrix that is compatible with other slice processing

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', true); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Binning
addOptional(p, 'binxy', 1);

% Start with registered or unregistered data
addOptional(p, 'useregistered', true);

% Which sections to use
addOptional(p, 'startsections', []);
addOptional(p, 'endsections', []);

% Ways to filter filter
addOptional(p, 'usefocusarea', false);
addOptional(p, 'usemedianthreshold', false);

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
% Get paths (run number does not matter)
spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);

% Load
load(fullfile(spcpaths.fp_out, spcpaths.mat), 'im2seg', 'im2seg2', 'icaguidata');

%% Preprocess ROIs
% Number of ROIs
nROI = length(icaguidata.CellAreas);

% Number of sections
nsections = spcpaths.n;

% Pass or fail
PassOrFail = ones(nROI, 1);

% Focus Areas
if p.usefocusarea
    % Get focus area
    FocusArea = getpoly(im2seg);
    
    % Out of focus area map
    OutOfFocusArea = icaguidata.AllFilters .* ~FocusArea;
    OutOfFocusAreaInds = unique(OutOfFocusArea(:));
    OutOfFocusAreaInds = OutOfFocusAreaInds(OutOfFocusAreaInds > 0);
    
    % Fail the areas that are out of focus
    PassOrFail(OutOfFocusAreaInds) = 0;
end

% Intensity threshold (median)
if p.usemedianthreshold
    % Find threshold
    int_thresh = median(im2seg2(:));
    
    % Loop through to check intensity
    for i = 1 : nROI
        if im2seg2(icaguidata.AllFilters == i) < int_thresh
            % Fail the cells that are dimmer than the threshold
            PassOrFail(i) = 0;
        end
    end
end

% Take a look at the new ROIs
if p.usefocusarea || p.usemedianthreshold
    % Make new ROIs
    newROI = zeros(size(icaguidata.AllFilters));
    
    for i = 1 : nROI
        if PassOrFail(i) > 0
            newROI(icaguidata.AllFilters == i) = i;
        end
    end
    
    % Take a look
    rgb = repmat(mat2gray(im2seg), [1 1 3]);
    rgb(:,:,3) = 0;
    rgb(:,:,2) = edge(newROI > 0);
    figure;imshow(rgb)
else
    newROI = icaguidata.AllFilters;
end

%% IO tiff files
% Load photon data
if p.useregistered
    % Load registered
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.regtif_photons));
else
    % Load unregistered
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_photons));
end

% Load tm data
if p.useregistered
    % Load registered
    im_tm = readtiff(fullfile(spcpaths.fp_out, spcpaths.regtif_tm));
else
    % Load unregistered
    im_tm = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_tm));
end

% Bin
if p.binxy > 1
    im_photon = binxy(im_photon, p.binxy);
    im_tm = binxy(im_tm, p.binxy);
end

%% Photons and Tm
% Initialze
Photons = zeros(nsections, nROI);
Tm = zeros(nsections, nROI);

% Dff values
Photonsdff = zeros(nsections, nROI);
Photonsdff_values = zeros(nROI, 1);
Photons_start = zeros(nROI, 1);
Tmdff_values = zeros(nROI, 1);

% Apply filters
for ROI_id = 1 : nROI
    if PassOrFail(ROI_id) > 0
        % Filter
        currfilt = icaguidata.AllFilters == ROI_id;

        % Area
        currarea = icaguidata.CellAreas(ROI_id);

        % Replicate filter for sections
        currfilt = repmat(currfilt, [1 1 nsections]);

        % Get photon counts
        trace_photons = squeeze(sum(sum(im_photon .* currfilt, 1), 2)) / currarea;
        Photons(:, ROI_id) = trace_photons;
            
        % Get Tm
        trace_tm = squeeze(sum(sum(im_tm .* currfilt, 1), 2)) / currarea;
        Tm(:, ROI_id) = trace_tm;
            
        % Photon dff
        Photonsdff(:, ROI_id) = trace_photons / mean(trace_photons(p.startsections));
        Photonsdff_values(ROI_id) = mean(trace_photons(p.endsections)) /...
            mean(trace_photons(p.startsections));
        Photons_start(ROI_id) = mean(trace_photons(p.startsections));
        
        % Tm dff
        Tmdff_values(ROI_id) = mean(trace_tm(p.endsections)) -...
            mean(trace_tm(p.startsections));
    end
end

%% Save to structure
% Make structures
OutputStruct.PassOrFail = PassOrFail;
OutputStruct.Areas = icaguidata.CellAreas;
OutputStruct.Photons = Photons;
OutputStruct.Tm = Tm;
OutputStruct.Photonsdff = Photonsdff;
OutputStruct.Photonsdff_values = Photonsdff_values;
OutputStruct.Photons_start = Photons_start;
OutputStruct.Tmdff_values = Tmdff_values;
OutputStruct.masks = icaguidata.AllFilters;
OutputStruct.masks_pass = newROI;

%% Save to file
save(fullfile(spcpaths.fp_out, spcpaths.mat), 'OutputStruct', '-append');
end