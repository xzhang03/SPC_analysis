function spcApplyROI(mouse, date, varargin)
% spcApplyROI applies ROIs to FLIM images

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Binning
addOptional(p, 'binxy', 1);

% Start with registered or unregistered data
addOptional(p, 'useregistered', true);

% Which sections to use (if numbers vary between experiments, use the
% highest number)
addOptional(p, 'nsections', 5);
addOptional(p, 'sections', []);

% Pick datasets
addOptional(p, 'dophotons', true); % Photon count
addOptional(p, 'dotm', true); % Mean lifetime
addOptional(p, 'plottm', true); % Make a plot in the end on tm data

% Save REF ROI
addOptional(p, 'saverefim', true);

% Save cross-run tm csv
addOptional(p, 'savexruntm', true);

% Multiple fovs (affects cross run names)
addOptional(p, 'multifov', false);
addOptional(p, 'fov', 1); % Specify fov

% Allowable pixels
addOptional(p, 'minallowabletm', 1000); % Minimal tm value that is allowed to be considered a cell (rejects 0 pixels)

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

% If not doing tm, not plotting tm
p.plottm = p.plottm & p.dotm;

%% IO
% Get paths (run number does not matter)
spcpaths = spcPath(mouse, date, 0, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit, 'multifov', p.multifov, 'fov', p.fov);

% Load
load(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_cell_clean', 'ROI_struct');

% Get the number of sections
nsections = p.nsections;

% Throw out the first section if
% p.sections is unspecified (First section contains frames where the
% resonant mirrows have yet to be up to speed)
if isempty(p.sections)
    p.sections = 2 : nsections;
end

% Get runs
runs = [ROI_struct.run];

% Get cells
ncells = size(ROI_cell_clean, 1);

%% Photons and Tm
% Areas
Area_mat = zeros(ncells, length(runs) + 1);

% Fill the first column with cell ID. Eventually it's cells x run. This
% matrix records all runs
Area_mat(:,1) = [ROI_cell_clean{:,1}];

% Photons initialize
if p.dophotons
    % Initialize
    photon_data_mat = zeros(ncells, nsections + 1);
    photon_avg_mat = zeros(ncells, length(runs) + 1);
        
    % Fill the first column with cell ID. Eventually it's cells x section.
    % This matrix changes for each run.
    photon_data_mat(:,1) = [ROI_cell_clean{:,1}];
    
    % Fill the first column with cell ID. Eventually it's cells x run. This
    % matrix records all runs
    photon_avg_mat(:,1) = [ROI_cell_clean{:,1}];
    
end

% Tm initialize
if p.dotm
    % Initialize
    tm_data_mat = zeros(ncells, nsections + 1);
    tm_avg_mat = zeros(ncells, length(runs) + 1);
    
    % Fill the first column with cell ID. Eventually it's cells x section .
    % This matrix changes for each run.
    tm_data_mat(:,1) = [ROI_cell_clean{:,1}];
    
    % Fill the first column with cell ID. Eventually it's cells x run. This
    % matrix records all runs
    tm_avg_mat(:,1) = [ROI_cell_clean{:,1}];
end


for run_ind = 1 : length(runs)
    % Get run path
    runpath = spcPath(mouse, date, runs(run_ind), 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Load photon data
    if p.dophotons && p.useregistered
        % Load registered
        [im_photon, nsections_photons] =...
            readtiff(fullfile(runpath.fp_out, runpath.regtif_photons));
    elseif p.dophotons
        % Load unregistered
        [im_photon, nsections_photons] =...
            readtiff(fullfile(runpath.fp_out, runpath.tif_photons));
    end
    
    % Load tm data
    if p.dotm && p.useregistered
        % Load registered
        [im_tm, nsections_tm] = readtiff(fullfile(runpath.fp_out, runpath.regtif_tm));
    elseif p.dotm
        % Load unregistered
        [im_tm, nsections_tm] = readtiff(fullfile(runpath.fp_out, runpath.tif_tm));
    end
    
    % Bin
    if p.dophotons && p.binxy > 1
        im_photon = binxy(im_photon, p.binxy);
    end
    if p.dotm && p.binxy > 1
        im_tm = binxy(im_tm, p.binxy);
    end
    
    % Apply shifts for xrun registration (does not exist for manual
    % matching)
    if ~isempty(ROI_struct(run_ind).shifts)
        if p.dophotons
            [~,im_photon] = stackRegisterMA_RR(im_photon, [], [], ...
                ones(size(im_photon,3),1) * ROI_struct(run_ind).shifts);
        end
        if p.dotm
            [~,im_tm] = stackRegisterMA_RR(im_tm, [], [], ...
                ones(size(im_tm,3),1) * ROI_struct(run_ind).shifts);
        end
    end
    
    % If doing tm, remove 0 pixels
    if p.dotm
        % Calculate an allowable region
        ROI_allow = min(im_tm, [], 3) >= p.minallowabletm;
    end
    
    % Apply filters
    for i = 1 : ncells
        % Filter id
        filterid = ROI_cell_clean{i,1};
        
        % Filter
        currfilt = ROI_struct(run_ind).ROI_xmatch_clean == filterid;
        
        % If doing tm, remove 0 pixels
        if p.dotm
            currfilt = currfilt .* ROI_allow;
        end
        
        % Area
        currarea = sum(currfilt(:));
        Area_mat(i, run_ind+1) = currarea;
        
        % Replicate filter for sections
        currfilt_photons = repmat(currfilt, [1 1 nsections_photons]);
        currfilt_tms = repmat(currfilt, [1 1 nsections_tm]);
        
        % Get photon counts
        if p.dophotons
            photon_data_mat(i, 2:end) =...
                vertcat(squeeze(sum(sum(im_photon .* currfilt_photons, 1), 2)) / currarea,...
                nan(nsections - nsections_photons, 1));
        end
        if p.dotm
            tm_data_mat(i, 2:end) =...
                vertcat(squeeze(sum(sum(im_tm .* currfilt_tms, 1), 2)) / currarea,...
                nan(nsections - nsections_tm, 1));
        end
    end
    
    % Find useful sections (not blurry or specified as not used here)
    sections_to_use = intersect(p.sections, ROI_struct(run_ind).sections);
    
    % Save data
    if p.dophotons
        % Update structure
        ROI_struct(run_ind).photon_stack = im_photon;
        ROI_struct(run_ind).photon_stack_med = median(im_photon(:,:,sections_to_use), 3);
        ROI_struct(run_ind).photon_mat = photon_data_mat;
        ROI_struct(run_ind).photon_vec = mean(photon_data_mat(:, sections_to_use + 1),2);
        
        % Collect all average data
        photon_avg_mat(:,run_ind + 1) = ROI_struct(run_ind).photon_vec;
    end
    if p.dotm
        % Update structure
        ROI_struct(run_ind).tm_stack = im_tm;
        ROI_struct(run_ind).tm_stack_med = median(im_tm(:,:,sections_to_use), 3);
        ROI_struct(run_ind).tm_mat = tm_data_mat;
        ROI_struct(run_ind).tm_vec = mean(tm_data_mat(:, sections_to_use + 1), 2);
        
        % Collect all average data
        tm_avg_mat(:,run_ind + 1) = ROI_struct(run_ind).tm_vec;
    end
end

%% Plot
if p.plottm
    figure
    hold on
    plot(1 : length(runs), tm_avg_mat(:,2:end), 'Color', [0.6 0.6 0.6]);
    plot(1 : length(runs), mean(tm_avg_mat(:,2:end),1), 'Color', [1 0 0], 'LineWidth', 5);
    hold off
    xlabel('Runs')
    ylabel('Tm (ps)')
end

%% Save to file
if ~p.dophotons
    photon_avg_mat = [];
end
if ~p.dotm
    tm_avg_mat = [];
end

% Save mat
save(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_struct', 'photon_avg_mat',...
    'tm_avg_mat', 'nsections', 'Area_mat', '-append')

if p.saverefim
    % Get centroids
    centroids = regionprops(ROI_struct(1).ROI_xmatch_clean, 'Centroid');
    cellids = cell2mat(ROI_cell_clean(:,1));
    centroids = centroids(cellids);
        
    % Make rgb
    rgb = repmat(mat2gray(ROI_struct(1).im),[1 1 3]);
    rgb(:,:,1) = (ROI_struct(1).ROI_xmatch_clean > 0) * 0.5;
    rgb(:,:,3) = 0;
    rgb = imresize(rgb, p.binxy);
    
    % show
    figure
    imshow(rgb)
    
    hold on
    % Loop and label
    for i = 1 : length(cellids)
        text(centroids(i).Centroid(1) * p.binxy, centroids(i).Centroid(2) * p.binxy,...
            num2str(cellids(i)), 'FontSize', 14, 'Color', [1 1 1]);
    end
    hold off
    
    % Save figure
    saveas(gcf, fullfile(spcpaths.fp_out, spcpaths.ROI_ref));
end

if p.savexruntm
    csvwrite(fullfile(spcpaths.fp_out, spcpaths.xruntm_csv), tm_avg_mat);
end

end