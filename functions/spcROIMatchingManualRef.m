function spcROIMatchingManualRef(mouse, date, runs, varargin)
% spcROImatchingRef matches ROIs between experiments, using another experiment
% as the main template. This function requires manual clicking.

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Figure variables
addOptional(p, 'figpos', []);
addOptional(p, 'redalpha', 0.5);

% Multiple fovs (affects cross run names)
addOptional(p, 'multifov', false);
addOptional(p, 'fov', 1); % Specify fov

% Reference info
addOptional(p, 'refmouse', '');
addOptional(p, 'refdate', '');
addOptional(p, 'refrun', 1);
addOptional(p, 'reffov', []);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Case and type
mouse = upper(mouse);

if isempty(p.refmouse)
    p.refmouse = mouse;
end

p.refmouse = upper(p.refmouse);

if ~ischar(date)
    date = num2str(date);
end

if ~ischar(p.refdate)
    p.refdate = num2str(p.refdate);
end

if isempty(p.reffov)
    p.reffov = p.fov;
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

% Figure locations
if isempty(p.figpos)
    figpos = [100 560 560 420; 100 50 560 420; 700 560 560 420; 700 50 560 420;...
        1300 560 560 420; 1300 50 560 420];
else
    figpos = p.figpos;
end

%% IO
% Initialize the ROIs (raw data, shifts, registered data, nICAs, areas)
ROI_struct = struct('run', [], 'sections', [], 'ROI_raw',[], 'im',[],...
    'ROI_xreg', [], 'ROI_xmatch', [], 'ROI_xmatch_clean', [], 'shifts', [],...
    'nROIs', NaN, 'Areas', []);
ROI_struct = repmat(ROI_struct, [length(runs), 1]);

% A cell to hold all the current ROIs
ROIs_left = cell(length(runs),1);

% A vector to hold the current number of ROIs
nROIs_left = zeros(length(runs),1);

% A cell to hold all the current RGBs (R = ROI, G = photon image)
RGB_cell = cell(length(runs),1);

for i = 1 : length(runs)
    % Get paths
    spcpaths = spcPath(mouse, date, runs(i), 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit, 'multifov', p.multifov, 'fov', p.fov);
    
    % Load
    loaded = load(fullfile(spcpaths.fp_out, spcpaths.mat), 'icaguidata', 'F', 'im2seg2');
    ROI_struct(i).run = runs(i);
    ROI_struct(i).sections = loaded.F;
    ROI_struct(i).ROI_raw = loaded.icaguidata.AllFilters; % ROIs
    ROI_struct(i).nROIs = length(loaded.icaguidata.CellAreas); % Number of ROIs
    ROI_struct(i).Areas = loaded.icaguidata.CellAreas; % Areas
    ROI_struct(i).im = loaded.im2seg2; % Original image
    
    ROI_struct(i).ROI_xmatch = zeros(size(loaded.icaguidata.AllFilters));   
    ROI_struct(i).ROI_xmatch_clean = ROI_struct(i).ROI_xmatch;
    
    % Current ROIs
    ROIs_left{i} = ROI_struct(i).ROI_raw;
    nROIs_left(i) = ROI_struct(i).nROIs;
    
    % RGB
    RGB = repmat(mat2gray(ROI_struct(i).im), [1 1 3]);
    RGB(:,:,1) = (ROI_struct(i).ROI_raw > 0) * p.redalpha;
    RGB(:,:,3) = 0;
    RGB_cell{i} = RGB;
end

%% IO reference
% Get paths
refpath = spcPath(p.refmouse, p.refdate, p.refrun, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit, 'multifov', p.multifov, 'fov', p.reffov);

% Indices for cells
cellref = load(fullfile(refpath.fp_out, refpath.xrun_mat), 'ROI_cell_clean');
cellref = cellref.ROI_cell_clean;
refinds = [cellref{:,1}];

% Refence image
imref = imread(fullfile(refpath.fp_out, refpath.ROI_ref));

% Get reference border
refbox = regionprops(imref(:,:,3) < 100, 'BoundingBox');
refbox = round(refbox.BoundingBox);

% Apply reference border
imref = imref(refbox(2) : refbox(2) + refbox(4), refbox(1) : refbox(1) + refbox(3), :); 

%% Intial plot
% Cells to get all the handles
hfigs = cell(length(runs),1);
hrgb = cell(length(runs),1);

% ROI cell
ROI_cell = cell(length(refinds), length(runs) + 1);
ROI_cell(:,1) = num2cell(refinds);

% Make reference figure
reffig = figure('name', sprintf('Reference'));
imshow(imref);
reftitle = title('Reference');
reftitle.Position(2) = 2;
set(reffig, 'Position', figpos(1,:));


for i = 1 : length(runs)
    % Figure
    hfigs{i} = figure('Position', figpos(i+1,:), 'name', sprintf('Run %i', runs(i)));
    hrgb{i} = imagesc(RGB_cell{i});
    title(sprintf('Run %i', runs(i)));
end

% Getting user input
for ii = 1 : length(refinds)
    % Reference index
    refind = refinds(ii);
    reftitle.String = sprintf('Find %i', refind);
    
    % User choice
    choice = questdlg(sprintf('Can you find %i?', refind), 'User input', 'Yes', 'No', 'Yes');
    
    switch choice
        case 'Yes'
            % Initialize chosen vector
            chosen = zeros(length(runs),1);
            
            % Get inputs
            for i = 1 : length(runs)
                % Get input
                figure(hfigs{i});
                
                % Figure which one was chosen
                while chosen(i) == 0
                    % Figure out which ROI
                    coord = round(ginput(1));
                    chosen(i) = ROIs_left{i}(coord(2), coord(1));
                end
                
                % Update cdata
                hrgb{i}.CData(:,:,1) = (ROIs_left{i} == chosen(i)) * p.redalpha;
            end
            
            % Apply inputs
            for i = 1 : length(runs)
                % Apply cross-matched ROIs
                ROI_struct(i).ROI_xmatch(ROIs_left{i} == chosen(i)) = refind;
                
                % Update cell
                ROI_cell{ii, i+1} = chosen(i);
            end
            
            % Remove ROI from pools
            for i = 1 : length(runs)
                % Decraese ROI count by 1
                nROIs_left(i) = nROIs_left(i) - 1;
                
                % Remove ROI
                ROIs_left{i}(ROIs_left{i} == chosen(i)) = 0;
                
                % Update RGB
                RGB_cell{i}(:,:,1) = (ROIs_left{i} > 0) * p.redalpha;
                
                % Update figure
                hrgb{i}.CData = RGB_cell{i};
            end
    end
end

%% Close
close(reffig);
for i = 1 : length(runs)
    close(hfigs{i});
end

%% Cleaning
% Matrix that denotes when ROIs are dropped
dropmat = cellfun(@isempty,ROI_cell);
stablecells = sum(dropmat,2) == 0;

% Clean ROI_cell
ROI_cell_clean = ROI_cell(stablecells,:);

% Clean ROIs
for i = 1 : size(ROI_cell_clean, 1)
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