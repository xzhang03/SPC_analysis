function spcROIMatchingManual2(mousecell, datecell, runs, varargin)
% spcROImatching matches ROIs between experiments, using the first
% experiment as the main template. This function requires manual clicking.
% This is a newer version that works with cellsort/cellpose

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Bin
addOptional(p, 'binxy', 1); % Bin used for generating the to-seg image before

% Figure variables
addOptional(p, 'figpos', []);
addOptional(p, 'redalpha', 0.5);

% String identifier (so other functions know which experiments are to be
% grouped together
addOptional(p, 'str', '');
addOptional(p, 'justchangestr', false);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Case and type
mousecell = upper(mousecell);

% User (add yourself if needed)
if isempty(p.user)
    switch mousecell(1:2)
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
    if length(mousecell) <= 6
        figpos = [100 560 560 420; 100 50 560 420; 700 560 560 420; 700 50 560 420;...
            1300 560 560 420; 1300 50 560 420];
    elseif length(mousecell) <= 12
        figpos = [100 750 400 300; 100 400 400 300; 100 50 400 300; 550 750 400 300;...
            550 400 400 300; 550 50 400 300; 1000 750 400 300; 1000 400 400 300;...
            1000 50 400 300; 1450 750 400 300; 1450 400 400 300; 1450 50 400 300];
    end
else
    figpos = p.figpos;
end

% Turn everything into cell
if ~iscell(mousecell)
    mousecell = {mousecell};
end

if ~iscell(datecell)
    datecell = {datecell};
end

if iscell(runs)
    runs = cell2mat(runs);
end
    
% Number of sets
nset = length(mousecell);

%% IO
% Initialize the ROIs (raw data, shifts, registered data, nICAs, areas)
ROI_struct = struct('mouse', [], 'date', [], 'run',[], ...
    'cellsort_spc', [], 'background', [], 'ROIs', [], 'ncells', []);
ROI_struct = repmat(ROI_struct, [nset, 1]);

% A cell to hold all the current ROIs
ROIs_left = cell(nset,1);

% A vector to hold the current number of ROIs
nROIs_left = zeros(nset,1);

% A cell to hold all the current RGBs (R = ROI, G = photon image)
RGB_cell = cell(nset,1);

hwait = waitbar(0, 'Loading');
for i = 1 : nset
    % waitbar
    waitbar(i/nset, hwait, sprintf('Loading %i/%i.', i, nset));
    
    % Basic values
    mouse = mousecell{i};
    date = datecell{i};
    runnum = runs(i);
    
    % Fill basic info
    ROI_struct(i).mouse = mouse;
    ROI_struct(i).date = date;
    ROI_struct(i).run = runnum;
    
    % Get paths
    spcpaths = spcPath(mouse, date, runnum, 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Load
    loaded = load(fullfile(spcpaths.fp_out, spcpaths.signals));
    
    % cell sort
    ROI_struct(i).cellsort_spc = loaded.cellsort_spc;
    ncells = length(loaded.cellsort_spc);
    ROI_struct(i).ncells = ncells;
    nROIs_left(i) = ncells;
    
    % Make a file of all ROIs
    ROIs = zeros(size(loaded.cellsort_spc(1).mask));
    if p.binxy > 1
        ROIs = binxy(ROIs, p.binxy);
    end
    
    for j = 1 : ncells
        mask = loaded.cellsort_spc(j).mask;
        if p.binxy > 1
            mask = binxy(mask, p.binxy);
        end
        ROIs(mask > 0) = j;
    end
    
    ROI_struct(i).ROIs = ROIs;
    ROIs_left{i} = ROIs;
    
    % Background
    background = loaded.movmean;
    if p.binxy > 1
        background = binxy(background, p.binxy);
    end
    ROI_struct(i).background = background;
        
    % RGB
    RGB = repmat(mat2gray(background), [1 1 3]);
    RGB(:,:,1) = (ROIs > 0) * p.redalpha;
    RGB(:,:,3) = 0;
    RGB_cell{i} = RGB;
end
close(hwait)

%% Intial plot
% Cells to get all the handles
hfigs = cell(length(runs),1);
hrgb = cell(length(runs),1);

% ROI cell
matchingmat = ones(ROI_struct(1).ncells, nset+1);
matchingmat(:,1) = 1 : ROI_struct(1).ncells;

for i = 1 : length(runs)
    % Figure
    hfigs{i} = figure('Position', figpos(i,:), 'name', sprintf('Run %i', runs(i)));
    hrgb{i} = imagesc(RGB_cell{i});
    title(sprintf('Date %s', datecell{i}));
    
    set(hfigs{i}, 'MenuBar', 'none');
    set(hfigs{i}, 'ToolBar', 'none');
end

% Current cell ID
cell_curr = 0;

% Getting user input
flagdone = false;

% Just changing strings
if p.justchangestr
    flagdone = true;
end

while ~flagdone
    % Initialize chosen vector
    chosen = zeros(nset,1);
    
    % If only 1 panel left, autochoose
    if sum(nROIs_left > 0) == 1
        panelleft = find(nROIs_left > 0);
        chosen(panelleft) = max(max(ROIs_left{panelleft}));
        autochoose = true;
    else
        autochoose = false;
    end
    
    % Advance cell count
    cell_curr = cell_curr + 1;
    
    % Get inputs
    for i = 1 : length(runs)
        % Get input
        figure(hfigs{i});
        
        % Skip
        skip = false;
        
        % Auto skip if none left
        if nROIs_left(i) == 0
            skip = true;
            chosen(i) = -1;
        end
        
        % Autochoose because only 1 panel left
        if autochoose
            skip = true;
        end
        
        % Figure which one was chosen
        while chosen(i) == 0
            % Figure out which ROI
            coord = round(ginput(1));
            
            if any(coord < 0) || coord(1)>hrgb{i}.XData(2) || coord(2)>hrgb{i}.YData(2)
                skip = true;
                chosen(i) = -1;
            else
                chosen(i) = ROIs_left{i}(coord(2), coord(1));
            end
        end

        % Update cdata
        if skip
             hrgb{i}.CData(:,:,1) = 0;
        else
             hrgb{i}.CData(:,:,1) = (ROIs_left{i} == chosen(i)) * p.redalpha;
        end
       
    end

    % Apply inputs
    for i = 1 : nset
        % Apply cellsort
        if chosen(i) > 0
            ROI_struct(i).cellsort_spc(chosen(i)).xrun_id = cell_curr;
        end
        
        % Update cell
        matchingmat(cell_curr, i+1)= chosen(i);
    end

    % Remove ROI from pools
    for i = 1 : length(runs)
        if chosen(i) > 0
            % Decraese ROI count by 1
            nROIs_left(i) = nROIs_left(i) - 1;

            % Remove ROI
            ROIs_left{i}(ROIs_left{i} == chosen(i)) = 0;

            % Update RGB
            RGB_cell{i}(:,:,1) = (ROIs_left{i} > 0) * p.redalpha;
        end
        
        % Update figure
        hrgb{i}.CData = RGB_cell{i};
    end

    % If no ROI left, it's done
    if sum(nROIs_left) == 0
        flagdone = true;
    end
end


%% Cleaning
% Matrix that denotes when ROIs are dropped
% Fix fitst column
if ~p.justchangestr
    matchingmat(:,1) = 1 : cell_curr;
end
for i = 1 : nset
    close(hfigs{i});
end

%% Saving
hwait = waitbar(0, 'Saving');
for i = 1 : nset
    % Save
    waitbar(i/nset, hwait, sprintf('Saving %i/%i.', i, nset));
    
    % Basic values
    mouse = mousecell{i};
    date = datecell{i};
    runnum = runs(i);
    
    % Get paths
    spcpaths = spcPath(mouse, date, runnum, 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Load
    loaded = load(fullfile(spcpaths.fp_out, spcpaths.signals));
    
    if ~p.justchangestr
        loaded.cellsort_spc = ROI_struct(i).cellsort_spc;
        loaded.xrun_ROIs = matchingmat;
    end
    loaded.fovstr = p.str;
    
    % Save
    save(fullfile(spcpaths.fp_out, spcpaths.signals), '-struct', 'loaded', '-v7.3');

end
disp('xrun matching done.')
close(hwait)

end