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

% Which sections to use (if not specified, it will be corrected to 2 : end
% below.
addOptional(p, 'sections', []);

% Pick datasets
addOptional(p, 'dophotons', true); % Photon count
addOptional(p, 'dotm', true); % Mean lifetime

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
spcpaths = spcPath(mouse, date, 0, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);

% Load
load(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_cell_clean', 'ROI_struct');

% Get the sections from the first experiment.
nsections = length(ROI_struct(1).sections);

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
        im_photon = readtiff(fullfile(runpath.fp_out, runpath.regtif_photons));
    elseif p.dophotons
        % Load unregistered
        im_photon = readtiff(fullfile(runpath.fp_out, runpath.tif_photons));
    end
    
    % Load tm data
    if p.dotm && p.useregistered
        % Load registered
        im_tm = readtiff(fullfile(runpath.fp_out, runpath.regtif_tm));
    elseif p.dotm
        % Load unregistered
        im_tm = readtiff(fullfile(runpath.fp_out, runpath.tif_tm));
    end
    
    % Bin
    if p.dophotons && p.binxy > 1
        im_photon = binxy(im_photon, p.binxy);
    end
    if p.dotm && p.binxy > 1
        im_tm = binxy(im_tm, p.binxy);
    end
    
    % Apply shifts for xrun registration
    if p.dophotons
        [~,im_photon] = stackRegisterMA_RR(im_photon, [], [], ...
            ones(size(im_photon,3),1) * ROI_struct(run_ind).shifts);
    end
    if p.dotm
        [~,im_tm] = stackRegisterMA_RR(im_tm, [], [], ...
            ones(size(im_tm,3),1) * ROI_struct(run_ind).shifts);
    end
    
    % Apply filters
    for i = 1 : ncells
        % Filter id
        filterid = ROI_cell_clean{i,1};
        
        % Filter
        currfilt = ROI_struct(run_ind).ROI_xmatch_clean == filterid;
        
        % Area
        currarea = sum(currfilt(:));
        
        % Replicate filter for sections
        currfilt = repmat(currfilt, [1 1 nsections]);
        
        % Get photon counts
        if p.dophotons
            photon_data_mat(i, 2:end) =...
                squeeze(sum(sum(im_photon .* currfilt, 1), 2)) / currarea;
        end
        if p.dotm
            tm_data_mat(i, 2:end) =...
                squeeze(sum(sum(im_tm .* currfilt, 1), 2)) / currarea;
        end
    end
    
    % Save data
    if p.dophotons
        % Update structure
        ROI_struct(run_ind).photon_stack = im_photon;
        ROI_struct(run_ind).photon_mat = photon_data_mat;
        ROI_struct(run_ind).photon_vec = mean(photon_data_mat(:, p.sections + 1),2);
        
        % Collect all average data
        photon_avg_mat(:,run_ind + 1) = ROI_struct(run_ind).photon_vec;
    end
    if p.dotm
        % Update structure
        ROI_struct(run_ind).tm_stack = im_tm;
        ROI_struct(run_ind).tm_mat = tm_data_mat;
        ROI_struct(run_ind).tm_vec = mean(tm_data_mat(:, p.sections + 1), 2);
        
        % Collect all average data
        tm_avg_mat(:,run_ind + 1) = ROI_struct(run_ind).tm_vec;
    end
end

%% Save to file
if ~p.dophotons
    photon_avg_mat = [];
end
if ~p.dotm
    tm_avg_mat = [];
end
save(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_struct', 'photon_avg_mat',...
    'tm_avg_mat', '-append')
end