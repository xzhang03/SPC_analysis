function spcROIMatching(mouse, date, runs, varargin)
% spcROImatching matches ROIs between experiments, using the first
% experiment as the main template

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Overlap threshold (In order for ROIs A and B to be considered the same,
% the overlapping area must be at least X% of A or X% of B
addOptional(p, 'OverlapThresh', 0.5);

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
% Initialize the ROIs (raw data, shifts, registered data, nICAs, areas)
ROI_struct = struct('run', [], 'sections', [], 'ROI_raw',[], 'ROI_xreg', [], ...
    'ROI_xmatch', [], 'ROI_xmatch_clean', [], 'shifts', [], 'nROIs', NaN, ...
    'Areas', []);
ROI_struct = repmat(ROI_struct, [length(runs), 1]);

for i = 1 : length(runs)
    % Get paths
    spcpaths = spcPath(mouse, date, runs(i), 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Load
    loaded = load(fullfile(spcpaths.fp_out, spcpaths.mat), 'icaguidata');
    ROI_struct(i).run = runs(i);
    ROI_struct(i).sections = spcpaths.cinds;
    ROI_struct(i).ROI_raw = loaded.icaguidata.AllFilters; % ROIs
    ROI_struct(i).nROIs = length(loaded.icaguidata.CellAreas); % Number of ROIs
    ROI_struct(i).Areas = loaded.icaguidata.CellAreas; % Areas
    ROI_struct(i).ROI_xmatch = zeros(size(loaded.icaguidata.AllFilters));   
    ROI_struct(i).ROI_xmatch_clean = ROI_struct(i).ROI_xmatch;
end

%% xy shift
% Loop through and shift
for i = 2 : length(runs)
    currentrun = runs(i);
    xyshift = stackRegisterMA_RR(ROI_struct(currentrun).ROI_raw, ROI_struct(1).ROI_raw);
    xyshift(:,3:4) = round(xyshift(:,3:4));
    [ROI_struct(i).shifts, ROI_struct(i).ROI_xreg] =...
        stackRegisterMA_RR(ROI_struct(currentrun).ROI_raw, [], [], xyshift);
    
    % Make sure it's rounded
    ROI_struct(i).ROI_xreg = round(ROI_struct(currentrun).ROI_xreg);
end

% No shift for reference
ROI_struct(1).shifts = [1 0 0 0];
ROI_struct(1).ROI_xreg = round(ROI_struct(1).ROI_raw);


%% ROI matching
% Initialize matrix
ROI_cell = cell(ROI_struct(1).nROIs, length(runs));
ROI_cell(:,1) = num2cell(1 : ROI_struct(1).nROIs);

for ROI_id = 1 : ROI_struct(1).nROIs
    % Reference ROI
    ref_roi = ROI_struct(3).ROI_xreg == ROI_id;
        
    for expt_id = 2 : length(runs)
        % Expt ROI
        expt_ROI = ROI_struct(expt_id).ROI_xreg .* ref_roi;
        
        % Overlap ROIs
        ovROIs = unique(expt_ROI(:));
        ovAreas = histc(expt_ROI(:), ovROIs);
        
        % Remove zeros
        ovAreas = ovAreas(ovROIs >= 1);
        ovROIs = ovROIs(ovROIs >= 1);
        
        % Area ratios
        ovAreaRatios_ref = ovAreas / ROI_struct(1).Areas(ROI_id);
        ovAreaRatios_expt = ovAreas ./ ROI_struct(expt_id).Areas(ovROIs);
        
        % Loop through
        for i = 1 : length(ovROIs)
            if ovAreaRatios_ref(i) >= p.OverlapThresh ||...
                    ovAreaRatios_expt(i) >= p.OverlapThresh
                % If the overlap is a significant part of the reference ROI
                % or of the experiment ROI
                % Label
                ROI_struct(expt_id).ROI_xmatch(ROI_struct(expt_id).ROI_xreg == ovROIs(i)) = ROI_id;
                
                % Update cell
                ROI_cell{ROI_id, expt_id} = [ROI_cell{ROI_id, expt_id}, ovROIs(i)];
            end
        end
    end
end

% No matching for reference
ROI_struct(1).ROI_xmatch = ROI_struct(1).ROI_xreg;


%% Cleaning
% Matrix that denotes when ROIs are dropped
dropmat = cellfun(@isempty,ROI_cell);
stablecells = sum(dropmat,2) == 0;

% Clean ROI_cell
ROI_cell_clean = ROI_cell(stablecells,:);

% Clean ROIs
for i = 1 : size(ROI_cell_clean)
    id = ROI_cell_clean{i,1};
    
    for expt_id = 1 : length(runs)
        % Label
        ROI_struct(expt_id).ROI_xmatch_clean(ROI_struct(expt_id).ROI_xmatch == id) = id;
    end
end

%% Saving
save(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_cell', 'ROI_struct',...
    'ROI_cell_clean', '-v7.3');
end