function spcROIMatchingManualSbx(mousecell_spc, datecell_spc, runs_spc,...
    mousecell_sbx, datecell_sbx, runs_sbx, optotunes_sbx, varargin)
% spcROImatching matches SPC ROIs onto sbx ROIs, using sbxROIs as templates.
% This function requires manual clicking.
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

% Sbx path varaibles
addOptional(p, 'useoptotune', false); % Whether using optotune or not
addOptional(p, 'pmt', 1);


% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Case and type
mousecell_spc = upper(mousecell_spc);

% User (add yourself if needed)
if isempty(p.user)
    switch mousecell_spc(1:2)
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
    if (length(mousecell_spc) + length(mousecell_sbx)) <= 6
        figpos = [100 560 560 420; 100 50 560 420; 700 560 560 420; 700 50 560 420;...
            1300 560 560 420; 1300 50 560 420];
    elseif (length(mousecell_spc) + length(mousecell_sbx)) <= 12
        figpos = [100 750 400 300; 100 400 400 300; 100 50 400 300; 550 750 400 300;...
            550 400 400 300; 550 50 400 300; 1000 750 400 300; 1000 400 400 300;...
            1000 50 400 300; 1450 750 400 300; 1450 400 400 300; 1450 50 400 300];
    end
else
    figpos = p.figpos;
end

% Turn everything into cell
if ~iscell(mousecell_spc)
    mousecell_spc = {mousecell_spc};
end

if ~iscell(datecell_spc)
    datecell_spc = {datecell_spc};
end

if iscell(runs_spc)
    runs_spc = cell2mat(runs_spc);
end
    
if iscell(runs_sbx)
    runs_sbx = cell2mat(runs_sbx);
end

% Number of sets
nset_spc = length(mousecell_spc);
    
% Number of sets
nset_sbx = length(mousecell_sbx);

% Just copying xrun id to xexpt id
if isempty(mousecell_spc)
    skipspc = true;
else 
    skipspc = false;
end
if isempty(mousecell_sbx)
    skipsbx = true;
else 
    skipsbx = false;
end

%% IO SPC
% Initialize the ROIs (raw data, shifts, registered data, nICAs, areas)
ROI_struct_spc = struct('mouse', [], 'date', [], 'run',[], ...
    'cellsort_spc', [], 'background', [], 'ROIs', [], 'ncells', []);
ROI_struct_spc = repmat(ROI_struct_spc, [nset_spc, 1]);

% A cell to hold all the current ROIs
ROIs_left_spc = cell(nset_spc,1);

% A vector to hold the current number of ROIs
nROIs_left_spc = zeros(nset_spc,1);

% A cell to hold all the current RGBs (R = ROI, G = photon image)
RGB_cell_spc = cell(nset_spc,1);

hwait = waitbar(0, 'Loading SPC');
for i = 1 : nset_spc
    % waitbar
    waitbar(i/nset_spc, hwait, sprintf('Loading SPC %i/%i.', i, nset_spc));
    
    % Basic values
    mouse = mousecell_spc{i};
    date = datecell_spc{i};
    runnum = runs_spc(i);
    
    % Fill basic info
    ROI_struct_spc(i).mouse = mouse;
    ROI_struct_spc(i).date = date;
    ROI_struct_spc(i).run = runnum;
    
    % Get paths
    spcpaths = spcPath(mouse, date, runnum, 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Load
    loaded_spc = load(fullfile(spcpaths.fp_out, spcpaths.signals));
    
    if i == 1
        matchingmat_spc = loaded_spc.xrun_ROIs;
    end
    
    % cell sort
    ROI_struct_spc(i).cellsort_spc = loaded_spc.cellsort_spc;
    ncells = length(loaded_spc.cellsort_spc);
    ROI_struct_spc(i).ncells = ncells;
    nROIs_left_spc(i) = ncells;
    
    % Make a file of all ROIs
    ROIs = zeros(size(loaded_spc.cellsort_spc(1).mask));
    if p.binxy > 1
        ROIs = binxy(ROIs, p.binxy);
    end
    
    for j = 1 : ncells
        mask = loaded_spc.cellsort_spc(j).mask;
        if p.binxy > 1
            mask = binxy(mask, p.binxy);
        end
        ROIs(mask > 0) = j;
    end
    
    ROI_struct_spc(i).ROIs = ROIs;
    ROIs_left_spc{i} = ROIs;
    
    % Background
    background = loaded_spc.movmean;
    if p.binxy > 1
        background = binxy(background, p.binxy);
    end
    ROI_struct_spc(i).background = background;
        
    % RGB
    RGB = repmat(mat2gray(background), [1 1 3]);
    RGB(:,:,1) = (ROIs > 0) * p.redalpha;
    RGB(:,:,3) = 0;
    RGB_cell_spc{i} = RGB;
end
close(hwait)

%% IO SBX
% Initialize
ROI_struct_sbx = struct('mouse', [], 'date', [], 'run',[], 'optotune',[],...
    'cellsort', [], 'background', [], 'ROIs', [], 'ncells', []);
ROI_struct_sbx = repmat(ROI_struct_sbx, [nset_sbx, 1]);

% A cell to hold all the current ROIs
ROIs_left_sbx = cell(length(runs_sbx),1);

% A vector to hold the current number of ROIs
nROIs_left_sbx = zeros(length(runs_sbx),1);

% A cell to hold all the current RGBs (R = ROI, G = photon image)
RGB_cell_sbx = cell(length(runs_sbx),1);

hwait = waitbar(0, 'Loading SBX');
for i = 1 : nset_sbx
    % waitbar
    waitbar(i/nset_sbx, hwait, sprintf('Loading SBX %i/%i.', i, nset_sbx));
    
    % Basic values
    mouse = mousecell_sbx{i};
    date = datecell_sbx{i};
    runnum = runs_sbx(i);
    optotune = optotunes_sbx(i);
    
    % Load
    ROI_struct_sbx(i).mouse = mouse;
    ROI_struct_sbx(i).date = date;
    ROI_struct_sbx(i).run = runnum;
    ROI_struct_sbx(i).optotune = optotune;
    
    % Signal file
    if p.useoptotune && ~isempty(optotune) && ~isnan(optotune)
        sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
    else
        sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
    end
    loaded_sbx = load(sigpath, '-mat', 'cellsort', 'xrun_ROIs');
    
    if i == 1
        matchingmat_sbx = loaded_sbx.xrun_ROIs;
    end
    
    % cell sort
    ROI_struct_sbx(i).cellsort = loaded_sbx.cellsort;
    ncells = length(loaded_sbx.cellsort);
    ROI_struct_sbx(i).ncells = ncells;
    nROIs_left_sbx(i) = ncells;
    
    % Make a file of all ROIs
    ROIs = zeros(size(loaded_sbx.cellsort(1).mask));
    for j = 1 : ncells
        ROIs(loaded_sbx.cellsort(j).mask) = j;
    end
    if p.binxy > 1
        ROIs = binxy(ROIs, p.binxy);
    end
    ROI_struct_sbx(i).ROIs = ROIs;
    ROIs_left_sbx{i} = ROIs;
    
    % background file
    [fp, ~, ~] = fileparts(sigpath);
    if ~p.useoptotune || isempty(optotune) || isnan(optotune)
        meanpath = fullfile(fp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, runnum));
    else
        meanpath = fullfile(fp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, runnum, optotune));
    end
    background = imread(meanpath);
    ROI_struct_sbx(i).background = background;
    
    % RGB
    RGB = repmat(mat2gray(background), [1 1 3]);
    RGB(:,:,1) = (ROIs > 0) * p.redalpha;
    RGB(:,:,3) = 0;
    RGB_cell_sbx{i} = RGB;
end
close(hwait);


%% Intial plot
% Cells to get all the handles
hfigs_spc = cell(length(runs_spc),1);
hrgb_spc = cell(length(runs_spc),1);

hfigs_sbx = cell(length(runs_sbx),1);
hrgb_sbx = cell(length(runs_sbx),1);

% ROI cell
if ~skipspc
    matchingmat_spc_xept = -ones(ROI_struct_spc(1).ncells, nset_spc+1);
    matchingmat_spc_xept(:,1) = 1 : ROI_struct_spc(1).ncells;
end
if ~skipsbx
    matchingmat_sbx_xept = -ones(ROI_struct_sbx(1).ncells, nset_sbx+1);
    matchingmat_sbx_xept(:,1) = 1 : ROI_struct_sbx(1).ncells;
end

for i = 1 : length(runs_spc)
    % Figure
    hfigs_spc{i} = figure('Position', figpos(i,:), 'name', sprintf('Run %i', runs_spc(i)));
    hrgb_spc{i} = imagesc(RGB_cell_spc{i});
    title(sprintf('SPC Date %s', datecell_spc{i}));
    
    set(hfigs_spc{i}, 'MenuBar', 'none');
    set(hfigs_spc{i}, 'ToolBar', 'none');
end

for i = 1 : length(runs_sbx)
    % Figure
    hfigs_sbx{i} = figure('Position', figpos(i + length(runs_spc),:), 'name',...
        sprintf('Run %i', runs_sbx(i)));
    hrgb_sbx{i} = imagesc(RGB_cell_sbx{i});
    title(sprintf('SBX Date %s', datecell_sbx{i}));
    
    set(hfigs_sbx{i}, 'MenuBar', 'none');
    set(hfigs_sbx{i}, 'ToolBar', 'none');
end


%% Clicking
% Current cell ID
cell_curr = 0;

% Getting user input
flagdone = false;

if skipspc || skipsbx
    flagdone = true;
end

while ~flagdone
    % Initialize chosen vector
    chosen_spc = zeros(nset_spc,1);
    chosen_sbx = zeros(nset_sbx,1);
%     % If only 1 panel left, autochoose
%     if sum(nROIs_left_spc > 0) == 1
%         panelleft = find(nROIs_left_spc > 0);
%         chosen(panelleft) = max(max(ROIs_left_spc{panelleft}));
%         autochoose = true;
%     else
%         autochoose = false;
%     end
    
    % Advance cell count
    cell_curr = cell_curr + 1;
    
    % Panel index spc
    if sum(nROIs_left_spc) > 0
        for i = 1 : nset_spc
            if all(chosen_spc <= 0)
                % Get input
                figure(hfigs_spc{i});

                % Skip
                skip_spc = false;

                % Auto skip if none left
                if nROIs_left_spc(i) == 0
                    skip_spc = true;
                    chosen_spc(i) = -1;
                end

        %         % Autochoose because only 1 panel left
        %         if autochoose
        %             skip = true;
        %         end

                % Figure which one was chosen
                while chosen_spc(i) == 0 && ~skip_spc
                    % Figure out which ROI
                    coord = round(ginput(1));
                    panelclicked = i;
                    
                    if any(coord < 0) || coord(1)>hrgb_spc{i}.XData(2) || coord(2)>hrgb_spc{i}.YData(2)
                        skip_spc = true;
                        chosen_spc(i) = -1;
                    else
                        chosen_spc(i) = ROIs_left_spc{i}(coord(2), coord(1));
                    end
                end
            end
        end
    end
    
    if any(chosen_spc > 0)
        xrun_id_spc = matchingmat_spc(:, panelclicked+1) == chosen_spc(chosen_spc > 0);
    else
        xrun_id_spc = [];
    end
    
    % Apply inputs
    for i = 1 : nset_spc
        % Current id
        curr_id = matchingmat_spc(xrun_id_spc, i+1);
        
        if curr_id > 0
            % Update UI
            hrgb_spc{i}.CData(:,:,1) = (ROIs_left_spc{i} == curr_id) * p.redalpha;
            
            % Update cellsort
            ROI_struct_spc(i).cellsort_spc(curr_id).xexpt_id = cell_curr;
            
            % Update cell
            matchingmat_spc_xept(cell_curr, i+1)= curr_id;
        else
            hrgb_spc{i}.CData(:,:,1) = 0;
        end
        
    end
    
    % Panel index sbx
    if sum(nROIs_left_sbx) > 0
        for i = 1 : nset_sbx
            if all(chosen_sbx <= 0)
                % Get input
                figure(hfigs_sbx{i});

                % Skip
                skip_sbx = false;

                % Auto skip if none left
                if nROIs_left_sbx(i) == 0
                    skip_sbx = true;
                    chosen_sbx(i) = -1;
                end

        %         % Autochoose because only 1 panel left
        %         if autochoose
        %             skip = true;
        %         end

                % Figure which one was chosen
                while chosen_sbx(i) == 0 && ~skip_sbx
                    % Figure out which ROI
                    coord = round(ginput(1));
                    panelclicked = i;

                    if any(coord < 0) || coord(1)>hrgb_sbx{i}.XData(2) || coord(2)>hrgb_sbx{i}.YData(2)
                        skip_sbx = true;
                        chosen_sbx(i) = -1;
                    else
                        chosen_sbx(i) = ROIs_left_sbx{i}(coord(2), coord(1));
                    end
                end
            end
        end
    end
    
    if any(chosen_sbx > 0)
        xrun_id_sbx = matchingmat_sbx(:, panelclicked+1) == chosen_sbx(chosen_sbx > 0);
    else
        xrun_id_sbx = [];
    end
    
    % Apply inputs
    for i = 1 : nset_sbx
        % Current id
        curr_id = matchingmat_sbx(xrun_id_sbx, i+1);
        
        if curr_id > 0
            % Update UI
            hrgb_sbx{i}.CData(:,:,1) = (ROIs_left_sbx{i} == curr_id) * p.redalpha;
            
            % Update cellsort
            ROI_struct_sbx(i).cellsort(curr_id).xexpt_id = cell_curr;
            
            % Update cell
            matchingmat_sbx_xept(cell_curr, i+1)= curr_id;
        else
            hrgb_sbx{i}.CData(:,:,1) = 0;
        end
        
    end

    % Remove ROI from spc pools
    for i = 1 : nset_spc
        % Current id
        curr_id = matchingmat_spc(xrun_id_spc, i+1);
        
        if curr_id > 0
            % Decraese ROI count by 1
            nROIs_left_spc(i) = nROIs_left_spc(i) - 1;

            % Remove ROI
            ROIs_left_spc{i}(ROIs_left_spc{i} == curr_id) = 0;

            % Update RGB
            RGB_cell_spc{i}(:,:,1) = (ROIs_left_spc{i} > 0) * p.redalpha;
        end
        
        % Update figure
        hrgb_spc{i}.CData = RGB_cell_spc{i};
    end
    
    % Remove ROI from sbx pools
    for i = 1 : nset_sbx
        % Current id
        curr_id = matchingmat_sbx(xrun_id_sbx, i+1);
        
        if curr_id > 0
            % Decraese ROI count by 1
            nROIs_left_sbx(i) = nROIs_left_sbx(i) - 1;

            % Remove ROI
            ROIs_left_sbx{i}(ROIs_left_sbx{i} == curr_id) = 0;

            % Update RGB
            RGB_cell_sbx{i}(:,:,1) = (ROIs_left_sbx{i} > 0) * p.redalpha;
        end
        
        % Update figure
        hrgb_sbx{i}.CData = RGB_cell_sbx{i};
    end

    % If no ROI left, it's done
    if sum(nROIs_left_spc) == 0 && sum(nROIs_left_sbx) == 0
        flagdone = true;
    end
end


%% Cleaning
% Matrix that denotes when ROIs are dropped
% Fix fitst column
if ~skipsbx
    matchingmat_sbx_xept(:,1) = 1 : size(matchingmat_sbx_xept,1);
end
if ~skipspc
    matchingmat_spc_xept(:,1) = 1 : size(matchingmat_spc_xept,1);
end

for i = 1 : nset_spc
    close(hfigs_spc{i});
end
for i = 1 : nset_sbx
    close(hfigs_sbx{i});
end

%% Copy xrun id to xexpt id
if skipspc
    % Copy sbx
    matchingmat_sbx_xept = matchingmat_sbx;
    
    for i = 1 : nset_sbx
        for j = 1:length(ROI_struct_sbx(i).cellsort)
            ROI_struct_sbx(i).cellsort(j).xexpt_id = ROI_struct_sbx(i).cellsort(j).xrun_id;
        end
    end
end

%% Saving
hwait = waitbar(0, 'Saving SPC');
for i = 1 : nset_spc
    % Save
    waitbar(i/nset_spc, hwait, sprintf('Saving SPC %i/%i.', i, nset_spc));
    
    % Basic values
    mouse = mousecell_spc{i};
    date = datecell_spc{i};
    runnum = runs_spc(i);
    
    % Get paths
    spcpaths = spcPath(mouse, date, runnum, 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Load
    loaded_spc = load(fullfile(spcpaths.fp_out, spcpaths.signals));
    
    loaded_spc.cellsort_spc = ROI_struct_spc(i).cellsort_spc;
    loaded_spc.xexpt_ROIs = matchingmat_spc_xept;
    
    % Save
    save(fullfile(spcpaths.fp_out, spcpaths.signals), '-struct', 'loaded_spc', '-v7.3');

end
close(hwait)

hwait = waitbar(0, 'Saving SBX');
for i = 1 : nset_sbx
    % Save
    waitbar(i/nset_sbx, hwait, sprintf('Saving SBX %i/%i.', i, nset_sbx));
    
    % Basic values
    mouse = mousecell_sbx{i};
    date = datecell_sbx{i};
    runnum = runs_sbx(i);
    optotune = optotunes_sbx(i);
    
    % Signal file
    if p.useoptotune && ~isempty(optotune) && ~isnan(optotune)
        sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
    else
        sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
    end
    
    % Update structures
    loaded_sbx = load(sigpath, '-mat');
    
    loaded_sbx.cellsort = ROI_struct_sbx(i).cellsort;
    loaded_sbx.xexpt_ROIs = matchingmat_sbx_xept;
    
    % Save
    save(sigpath, '-struct', 'loaded_sbx', '-v7.3');
end
close(hwait)

disp('xexpt matching done.')
end