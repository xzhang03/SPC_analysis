function spcNeuropil(mouse, date, varargin)
% spcNeurpil calculate neuropil lifetime as the average or weighted average
% of lifetimes in the neuropil ring defined by a inner distance and an outter
% distance.

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Binning
addOptional(p, 'binxy', 1);

% Which sections to use
addOptional(p, 'sections', []);

% Pick datasets
addOptional(p, 'dophotons', true); % Photon count
addOptional(p, 'dotm', true); % Mean lifetime
addOptional(p, 'plottm', true); % Make a plot in the end on tm data

% Ring
addOptional(p, 'innerspacing', 2); % pixels after binning
addOptional(p, 'outterspacing', 5);% pixels after binnng
addOptional(p, 'minpixels', 20); % Minimal number of pixels for neuropil calculations (this number is after binning)
addOptional(p, 'excludeallROIs', true); % Exclude all ROIs or only the cleaned ROIs

% Calculation method
addOptional(p, 'averaging_method', 'weighted'); % Can be 'weighted' or 'uniform'

% Save REF ROI
addOptional(p, 'saverefim', true);

% Save cross-run tm csv
addOptional(p, 'savexruntm', true);

% Multiple fovs (affects cross run names)
addOptional(p, 'multifov', false);
addOptional(p, 'fov', 1); % Specify fov

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

% Cannot do tm without photondata (since tm is weighted average by photon
% count)
p.dophotons = p.dophotons | p.dotm;

%% IO
% Get paths (run number does not matter)
spcpaths = spcPath(mouse, date, 0, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit, 'multifov', p.multifov, 'fov', p.fov);

% Load
load(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_cell_clean', 'ROI_struct', 'nsections');

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
Areanp_mat = zeros(ncells, length(runs) + 1);

% Fill the first column with cell ID. Eventually it's cells x run. This
% matrix records all runs
Areanp_mat(:,1) = [ROI_cell_clean{:,1}];

if p.dophotons
    % Initialize
    photon_datanp_mat = zeros(ncells, nsections + 1);
    photon_avgnp_mat = zeros(ncells, length(runs) + 1);
        
    % Fill the first column with cell ID. Eventually it's cells x section.
    % This matrix changes for each run.
    photon_datanp_mat(:,1) = [ROI_cell_clean{:,1}];
    
    % Fill the first column with cell ID. Eventually it's cells x run. This
    % matrix records all runs
    photon_avgnp_mat(:,1) = [ROI_cell_clean{:,1}];
    
end
if p.dotm
    % Initialize
    tm_datanp_mat = zeros(ncells, nsections + 1);
    tm_avgnp_mat = zeros(ncells, length(runs) + 1);
    
    % Fill the first column with cell ID. Eventually it's cells x section .
    % This matrix changes for each run.
    tm_datanp_mat(:,1) = [ROI_cell_clean{:,1}];
    
    % Fill the first column with cell ID. Eventually it's cells x run. This
    % matrix records all runs
    tm_avgnp_mat(:,1) = [ROI_cell_clean{:,1}];
end

for run_ind = 1 : length(runs)
    % Get the images
    if p.dophotons
        % Stack
        im_photon  = ROI_struct(run_ind).photon_stack;
        nsections_photons = size(im_photon, 3);
        
        % Convert to 2D if using weighted method
        if strcmp(p.averaging_method, 'weighted')
            im_photon = ...
                reshape(im_photon, [size(im_photon,1) * size(im_photon,2), nsections_photons]);
        end
    end
    if p.dotm
        % Stack
        im_tm  = ROI_struct(run_ind).tm_stack;
        nsections_tm = size(im_tm, 3);
        
        % Convert to 2D if using weighted method
        if strcmp(p.averaging_method, 'weighted')
            im_tm = ...
                reshape(im_tm, [size(im_tm,1) * size(im_tm,2), nsections_tm]);
        end
    end
    
    % Initialize np image
    ROInp_xmatch_clean = zeros(size(ROI_struct(run_ind).ROI_xmatch_clean));
    
    % A negative image to subtract from neuropil
    if p.excludeallROIs
        % Subtract all ROIs
        np_negative = ROI_struct(run_ind).ROI_raw > 0;
    else
        % Subtract only the clean ROIs
        np_negative = ROI_struct(run_ind).ROI_xmatch_clean > 0;
    end
    
    % Loop through cells
    for cellind = 1 : size(ROI_cell_clean,1)
        % Cell id
        id = ROI_cell_clean{cellind, 1};
        
        % Make neuropil
        [np, np_size] = MakeNPRing(ROI_struct(run_ind).ROI_xmatch_clean == id,...
            p.innerspacing, p.outterspacing, np_negative, p.minpixels);
        Areanp_mat(cellind, run_ind + 1) = np_size;
        
        % Log neuropil
        ROInp_xmatch_clean = ROInp_xmatch_clean + np;
        
        switch p.averaging_method
            case 'uniform'
                % Replicate filter for sections
                np_photons = repmat(np, [1 1 nsections_photons]);
                np_tms = repmat(np, [1 1 nsections_tm]);

                % Get photon counts
                if p.dophotons
                    photon_datanp_mat(cellind, 2:end) =...
                        vertcat(squeeze(sum(sum(im_photon .* np_photons, 1), 2)) / np_size,...
                        nan(nsections - nsections_photons, 1));
                end
                if p.dotm
                    tm_datanp_mat(cellind, 2:end) =...
                        vertcat(squeeze(sum(sum(im_tm .* np_tms, 1), 2)) / np_size,...
                        nan(nsections - nsections_tm, 1));
                end
            case 'weighted'
                % Replicate filter for sections
                np_photons = (np(:) * ones(1, nsections_photons)) > 0;
                np_tm = (np(:) * ones(1, nsections_tm)) > 0;
                
                % Photons
                if p.dophotons
                    % Apply filters
                    mat_photons = reshape(im_photon(np_photons), [np_size, nsections_photons]);
                    photon_datanp_mat(cellind, 2:end) = ...
                        horzcat(mean(mat_photons,1), nan(1, nsections - nsections_photons));
                end
                
                % Photons
                if p.dotm
                    % Apply filters
                    mat_tm = reshape(im_tm(np_tm), [np_size, nsections_tm]);
                    tm_datanp_mat(cellind, 2:end) = ...
                        horzcat(sum(mat_photons .* mat_tm,1) ./ sum(mat_photons, 1),...
                        nan(1, nsections - nsections_photons));
                end
             
        end
    end
    
    % Find useful sections (not blurry or specified as not used here)
    sections_to_use = intersect(p.sections, ROI_struct(run_ind).sections);
    
    % Save data
    if p.dophotons
        % Update structure
        ROI_struct(run_ind).photon_np_mat = photon_datanp_mat;
        ROI_struct(run_ind).photon_np_vec = mean(photon_datanp_mat(:, sections_to_use + 1),2);
        
        % Collect all average data
        photon_avgnp_mat(:,run_ind + 1) = ROI_struct(run_ind).photon_np_vec;
    end
    if p.dotm
        % Update structure
        ROI_struct(run_ind).tm_np_mat = tm_datanp_mat;
        ROI_struct(run_ind).tm_np_vec = mean(tm_datanp_mat(:, sections_to_use + 1), 2);
        
        % Collect all average data
        tm_avgnp_mat(:,run_ind + 1) = ROI_struct(run_ind).tm_np_vec;
    end
    
    % Save neuropil ROIs
    ROI_struct(run_ind).ROInp_xmatch_clean = ROInp_xmatch_clean;
end

%% Plot
if p.plottm
    figure
    hold on
    plot(1 : length(runs), tm_avgnp_mat(:,2:end), 'Color', [0.6 0.6 0.6]);
    plot(1 : length(runs), mean(tm_avgnp_mat(:,2:end),1), 'Color', [1 0 0], 'LineWidth', 5);
    hold off
    xlabel('Runs')
    ylabel('Neuropil Tm (ps)')
end

%% Save to file
if ~p.dophotons
    photon_avgnp_mat = [];
end
if ~p.dotm
    tm_avgnp_mat = [];
end

% Save mat
save(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_struct', 'photon_avgnp_mat',...
    'tm_avgnp_mat', '-append')

if p.saverefim
    % Get centroids
    centroids = regionprops(ROI_struct(1).ROI_xmatch_clean, 'Centroid');
    cellids = cell2mat(ROI_cell_clean(:,1));
    centroids = centroids(cellids);
        
    % Make rgb (with the first run)
    rgb = repmat(mat2gray(ROI_struct(1).im),[1 1 3]);
    rgb(:,:,1) = ROI_struct(1).ROInp_xmatch_clean * 0.5;
    rgb(:,:,3) = (ROI_struct(1).ROI_xmatch_clean > 0) * 0.5;
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
    saveas(gcf, fullfile(spcpaths.fp_out, spcpaths.ROInp_ref));
end

if p.savexruntm
    csvwrite(fullfile(spcpaths.fp_out, spcpaths.xruntmnp_csv), tm_avgnp_mat);
end
end