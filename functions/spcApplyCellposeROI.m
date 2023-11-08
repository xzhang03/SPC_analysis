function spcApplyCellposeROI(mouse, date, run, varargin)
% spcApplySbxROI applies Sbx ROIs to FLIM images

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Force
addOptional(p, 'force', false); % Force

% Registration variables
addOptional(p, 'sourcetype', 'demonsreg'); % Input type can be 'warped' or 'demonsreg'

% Binning
addOptional(p, 'binxy', 1);

% Tm
addOptional(p, 'mintm', 1000); % Min tm

% IEM
addOptional(p, 'miniem', 300); % Min tm

% Threshold
addOptional(p, 'minarea', 20); % Min number of pixels for a cell (after binning)

% Neuropil
addOptional(p, 'npsize', [14, 6]); % Neuropil size, avoidance size
addOptional(p, 'minnparea', 100); % Minimal number of pixels for neuropil

% Pick datasets
addOptional(p, 'dophotons', true); % Photon count
addOptional(p, 'dotm', true); % Mean lifetime
addOptional(p, 'doiem', true); % IEM lifetime
addOptional(p, 'plotscatter', true); % Make a plot in the end on tm data

% Frame info for dff and dt
addOptional(p, 'StartFrame', 1); % Which frame to start
addOptional(p, 'StartFrameN', 12); % How many frames to read for the initial stack
addOptional(p, 'EndFrame', []); % Which frame to begin the end stack. Leave empty to read the last X frames
addOptional(p, 'EndFrameN', 12); % How many frames to read for the end stack

% Save REF ROI
addOptional(p, 'saverefim', true);

% Save cross-run tm csv
addOptional(p, 'savetmcsv', true);

% GRIN
addOptional(p, 'GRIN', false); % Remove anything that is not within the GRIN surface
addOptional(p, 'reusegrinface', true); % Try to reuse grin face from a previous run
addOptional(p, 'reusegrinfacerun', []); % The rune to use grin face from
addOptional(p, 'grinface', []); % Pass the face of the GRIN lens here if you don't want manual selection

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


%% IO paths
% Get flim paths
spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);

% Image paths
% Switch types
switch p.sourcetype
    case 'raw'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.tif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.tif_tm);
        fp_iem = fullfile(spcpaths.fp_out, spcpaths.tif_iem);
    case 'registered'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.regtif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.regtif_tm);
        fp_iem = fullfile(spcpaths.fp_out, spcpaths.regtif_iem);
    case 'warped'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.warptif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.warptif_tm);
        fp_iem = fullfile(spcpaths.fp_out, spcpaths.warptif_iem);
    case 'demonsreg'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.demregtif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.demregtif_tm);
        fp_iem = fullfile(spcpaths.fp_out, spcpaths.demregtif_iem);
end


%% Load data
% Check exist
sigpath = fullfile(spcpaths.fp_out, spcpaths.signals);
if exist(sigpath, 'file') && ~p.force
    redo = input('Signal file already exist, redo? (1 = yes, 0 - no): ');
    if redo ~= 1
        return;
    end
end

% Load the masks
maskpath = fullfile(spcpaths.fp_out, spcpaths.cp_masks);
masks = imread(maskpath);

% Initialize a place to hold all the neuropils
npall = zeros(size(masks));

% Load background
tosegpath = fullfile(spcpaths.fp_out, spcpaths.cp_toseg);
movmean = readtiff(tosegpath);

% Basic paramteres
ncells = double(max(masks(:)));

% Get the face of GRIN
if p.GRIN
    if isempty(p.grinface)
        if p.reusegrinface 
            if isempty(p.reusegrinfacerun) && exist(sigpath, 'file')
                % Reuse grin face from a previous analysis
                loaded = load(sigpath, '-mat', 'psig');
                p.grinface = loaded.psig.grinface;
            elseif ~isempty(p.reusegrinfacerun)
                % Reuse grin face from a different run
                spcpath_grin = spcPath(mouse, date, p.reusegrinfacerun, 'server', p.server, 'user', p.user,...
                    'slice', p.slice, 'cdigit', p.cdigit);
                loaded = load(fullfile(spcpath_grin.fp_out, spcpath_grin.signals), '-mat', 'psig');
                p.grinface = loaded.psig.grinface;
            else
                p.grinface = getpoly(movmean, 'Select the face of the GRIN lens.');
                p.grinface = p.grinface;
            end
        else
            p.grinface = getpoly(movmean, 'Select the face of the GRIN lens.');
            p.grinface = p.grinface;
        end
    end
    
    % Update masks
    masks = masks .* uint16(p.grinface);
end



%% Initialize
% datastruct
trace_ini = struct('raw', [], 'neuropil', [], 'subtracted', [], 'raw_delta', [],...
    'neuropil_delta', [], 'subtracted_delta', []); % Delta is dff or dt
cellsort_spc = struct('id', [], 'mask', [], 'area', [], 'neuropil', [], 'neuropilarea', [], 'group_number', [],...
    'photon_trace', trace_ini, 'tm_trace', trace_ini, 'tm_start',[], 'iem_trace', trace_ini, 'iem_start',[]);
cellsort_spc = repmat(cellsort_spc, [ncells, 1]);

%% Photons and Tm
% Load photon data
if p.dophotons && ~exist('im_photon', 'var')
    % Load warped
    [im_photon, ~] = readtiff(fp_photon);
end

% Load tm data
if p.dotm
    if exist(fp_tm, 'file')
        % Load warped
        [im_tm, ~] = readtiff(fp_tm);   
    else
        disp('No Tm source file. Skip.');
        p.dotm = false;
    end
end

% Load IEM data
if p.doiem
    if exist(fp_iem, 'file')
        % Load warped
        [im_iem, ~] = readtiff(fp_iem);   
    else
        disp('No IEM source file. Skip.');
        p.doiem = false;
    end
end

% Get the number of sections
if exist('im_photon', 'var')
    nsections = size(im_photon, 3);
elseif exist('im_tm', 'var')
    nsections = size(im_tm, 3);
end

% Fix endframe
if isempty(p.EndFrame)
    p.EndFrame = nsections - p.EndFrameN + 1;
end

% Bin
if p.dophotons && p.binxy > 1
    im_photon = binxy(im_photon, p.binxy);
end
if p.dotm && p.binxy > 1
    im_tm = binxy(im_tm, p.binxy);
end
if p.doiem && p.binxy > 1
    im_iem = binxy(im_iem, p.binxy);
end

if p.dotm
    % Get rid of bad pixels (tm way too low)
    goodpixels_tm = im_tm >= p.mintm;
    % goodpixels = mean(goodpixels,3);
    im_tm(~goodpixels_tm) = nan;
end
if p.doiem
    % Get rid of bad pixels (tm way too low)
    goodpixels_iem = im_iem >= p.miniem;
    % goodpixels = mean(goodpixels,3);
    im_iem(~goodpixels_iem) = nan;
end

% Neuropil elements
npelbig = strel('disk', p.npsize(1));
npelsmall = strel('disk', p.npsize(2));

% Index for real cells
irealcell = 0;

% All ROIs
allmasks = imdilate(masks, npelsmall) > 0;

if p.dotm || p.dophotons
    % Apply filters
    hwait = waitbar(0, 'Applying masks');
    for i = 1 : ncells
        if mod(i,5) == 0
            waitbar(i/ncells, hwait, sprintf('Applying masks %i/%i', i, ncells));
        end

        % Current mask
        maskcurr = masks == i;
        maskarea = sum(maskcurr(:));
        
        % Check area
        if maskarea >= p.minarea
            % Advance real cell counter
            irealcell = irealcell + 1;
            
            % Neuropil (before binning)
            npcurr = imdilate(maskcurr, npelbig);
            if p.GRIN
                npcurr = npcurr .* double(p.grinface);
            end
            npcurr = npcurr .* ~allmasks;
            nparea = sum(npcurr(:) > 0);
            
            % Make sure np areai is ok
            if nparea < p.minnparea
                npsizetry = p.npsize(1);
                
                while nparea < p.minnparea
                    npsizetry = npsizetry + 1;
                    npcurr = imdilate(maskcurr, strel('disk', npsizetry));
                    
                    if p.GRIN
                        npcurr = npcurr .* p.grinface;
                    end
                    npcurr = npcurr .* ~allmasks;
                    nparea = sum(npcurr(:) > 0);
                end
            end
            
            % Save id
            cellsort_spc(irealcell).id = i;
            
            % Save mask
            cellsort_spc(irealcell).mask = maskcurr;

            % Save area
            cellsort_spc(irealcell).area = maskarea;
            
            % Save neuropil
            cellsort_spc(irealcell).neuropil = npcurr;
            npall = npall + npcurr;
            
            % Save neuropil area
            cellsort_spc(irealcell).neuropilarea = nparea;

            % Bin filter if needed
            if p.binxy > 1
                maskcurr = binxy(maskcurr, p.binxy);
                npcurr = binxy(npcurr, p.binxy);
                maskarea = maskarea / p.binxy^2;
                nparea = nparea / p.binxy^2;
            end
            
            % Replicate filter for sections
            mask_stack = repmat(maskcurr, [1 1 nsections]) > 0;
            np_stack = repmat(npcurr, [1 1 nsections]) > 0;
            
            % Get photon counts
            if p.dophotons
                % raw photon trace
                photon_trace =...
                    squeeze(nansum(nansum(im_photon .* mask_stack, 1), 2)) / maskarea;
                cellsort_spc(irealcell).photon_trace.raw = photon_trace;
                
                % np photon trace
                photon_np_trace = ...
                    squeeze(nansum(nansum(im_photon .* np_stack, 1), 2)) / nparea;
                cellsort_spc(irealcell).photon_trace.neuropil = photon_np_trace;
                
                % subtracted photon trace
                photon_subtracted_trace = photon_trace - photon_np_trace;
                cellsort_spc(irealcell).photon_trace.subtracted = photon_subtracted_trace;
                
                % dff photon trace
                dff_trace =...
                    photon_trace / mean(photon_trace(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
                cellsort_spc(irealcell).photon_trace.raw_delta = dff_trace;
                
                % dff photon np trace
                dff_np_trace =...
                    photon_np_trace / mean(photon_np_trace(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
                cellsort_spc(irealcell).photon_trace.neuropil_delta = dff_np_trace;
                
                % dff photon subtracted trace
                dff_subtracted_trace =...
                    photon_subtracted_trace / mean(photon_subtracted_trace(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
                cellsort_spc(irealcell).photon_trace.subtracted_delta = dff_subtracted_trace;
            end
            
            if p.dotm
                % raw tm trace
                tm_trace =...
                    squeeze(nansum(nansum(im_tm .* mask_stack, 1), 2));
                tm_trace_area = squeeze(nansum(nansum(goodpixels_tm .* mask_stack, 1), 2));
                tm_trace = tm_trace ./ tm_trace_area;
                cellsort_spc(irealcell).tm_trace.raw = tm_trace;
                
                % np tm trace
                tm_np_trace = ...
                    squeeze(nansum(nansum(im_tm .* np_stack, 1), 2));
                tm_np_trace_area = squeeze(nansum(nansum(goodpixels_tm .* np_stack, 1), 2));
                tm_np_trace = tm_np_trace ./ tm_np_trace_area;
                cellsort_spc(irealcell).tm_trace.neuropil = tm_np_trace;
                
                % subtracted tm trace
                tm_subtracted_trace = tm_trace - tm_np_trace;
                cellsort_spc(irealcell).tm_trace.subtracted = tm_subtracted_trace;

                % dt tm trace
                dt_trace =...
                    tm_trace - mean(tm_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
                cellsort_spc(irealcell).tm_trace.raw_delta = dt_trace;
                
                % dt tm np trace
                dt_np_trace =...
                    tm_np_trace - mean(tm_np_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
                cellsort_spc(irealcell).tm_trace.neuropil_delta = dt_np_trace;
                
                % dt tm subtracted trace
                dt_subtracted_trace =...
                    tm_subtracted_trace - mean(tm_subtracted_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
                cellsort_spc(irealcell).tm_trace.subtracted_delta = dt_subtracted_trace;
                
                % Tm start
                cellsort_spc(irealcell).tm_start = mean(tm_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
            end
            
            if p.doiem
                % raw iem trace
                iem_trace =...
                    squeeze(nansum(nansum(im_iem .* mask_stack, 1), 2));
                iem_trace_area = squeeze(nansum(nansum(goodpixels_iem .* mask_stack, 1), 2));
                iem_trace = iem_trace ./ iem_trace_area;
                cellsort_spc(irealcell).iem_trace.raw = iem_trace;
                
                % np iem trace
                iem_np_trace = ...
                    squeeze(nansum(nansum(im_iem .* np_stack, 1), 2));
                iem_np_trace_area = squeeze(nansum(nansum(goodpixels_iem .* np_stack, 1), 2));
                iem_np_trace = iem_np_trace ./ iem_np_trace_area;
                cellsort_spc(irealcell).iem_trace.neuropil = iem_np_trace;
                
                % subtracted iem trace
                iem_subtracted_trace = iem_trace - iem_np_trace;
                cellsort_spc(irealcell).iem_trace.subtracted = iem_subtracted_trace;

                % dt iem trace
                dt_trace =...
                    iem_trace - mean(iem_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
                cellsort_spc(irealcell).iem_trace.raw_delta = dt_trace;
                
                % dt iem np trace
                dt_np_trace =...
                    iem_np_trace - mean(iem_np_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
                cellsort_spc(irealcell).iem_trace.neuropil_delta = dt_np_trace;
                
                % dt iem subtracted trace
                dt_subtracted_trace =...
                    iem_subtracted_trace - mean(iem_subtracted_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
                cellsort_spc(irealcell).iem_trace.subtracted_delta = dt_subtracted_trace;
                
                % IEM start
                cellsort_spc(irealcell).iem_start = mean(iem_trace(p.StartFrame : p.StartFrame+p.StartFrameN-1));
            end
        end
    end
    close(hwait)
end

% Shrink cellsort
nrealcell = irealcell;
cellsort_spc = cellsort_spc(1:nrealcell);

%% Plot
if p.plotscatter
    % Initialize
    photonsvec = nan(nrealcell,1);
    if p.dotm
        tmvec = nan(nrealcell,1);
    end
    if p.doiem
        iemvec = nan(nrealcell,1);
    end
    
    for i = 1 : nrealcell
        % Load
        photonsvec(i) =...
            mean(cellsort_spc(i).photon_trace.subtracted_delta(p.EndFrame:p.EndFrame+p.EndFrameN-1))...
            - mean(cellsort_spc(i).photon_trace.subtracted_delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        
        if p.dotm
            tmvec(i) =...
                mean(cellsort_spc(i).tm_trace.subtracted_delta(p.EndFrame:p.EndFrame+p.EndFrameN-1))...
                - mean(cellsort_spc(i).tm_trace.subtracted_delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        end
        
        if p.doiem
            iemvec(i) =...
                mean(cellsort_spc(i).iem_trace.subtracted_delta(p.EndFrame:p.EndFrame+p.EndFrameN-1))...
                - mean(cellsort_spc(i).iem_trace.subtracted_delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        end
    end

    if p.dotm
        figure
        scatter(photonsvec, tmvec);
        title('Tm vs photons')
    end
    if p.doiem
        figure
        scatter(photonsvec, iemvec);
        title('IEM vs photons')
    end
end

%% Save to file
% Save structure
sigstruct = struct('cellsort_spc', cellsort_spc, 'psig', p, 'movmean', movmean,...
    'masks', masks, 'neuropil', npall);

% Save file
save(sigpath, '-struct', 'sigstruct', '-v7.3');

if p.saverefim
    % Make rgb
    maskrgb = repmat(mat2gray(movmean), [1 1 3]);
    ROIs_bin = mat2gray(masks >= 1);
    maskrgb(:,:,1) = maskrgb(:,:,1) + ROIs_bin;
    
    % Show rgb
    figure
    imshow(maskrgb);
    
    % Label
    for i = 1 : ncells
        if sum(sum(masks == i)) > 0
            cen = regionprops(masks == i, 'Centroid');
            cen = cen.Centroid;
            text(cen(1), cen(2), num2str(i), 'Color', [0 0 0], 'FontSize', 14);
        end
    end
    
    saveas(gcf, fullfile(spcpaths.fp_out, spcpaths.run_ROI_ref));
    
    % Make neuropil rgb
    nprgb = repmat(mat2gray(movmean), [1 1 3]);
    ROIs_bin = mat2gray(npall >= 1);
    nprgb(:,:,1) = nprgb(:,:,1) + ROIs_bin;

    % Show rgb
    figure
    imshow(nprgb);
    
    % Label
    for i = 1 : ncells
        if sum(sum(masks == i)) > 0
            cen = regionprops(masks == i, 'Centroid');
            cen = cen.Centroid;
            text(cen(1), cen(2), num2str(i), 'Color', [0 0 0], 'FontSize', 14);
        end
    end
    
    saveas(gcf, fullfile(spcpaths.fp_out, spcpaths.run_npROI_ref));
end

if p.savetmcsv
    % Iinitialize
    tmmat = nan(nrealcell,5);
    
    for i = 1 : nrealcell
        % Tm
        if p.dotm
            % Load
            tmmat(i,3) = mean(cellsort_spc(i).tm_trace.subtracted_delta(p.EndFrame:p.EndFrame+p.EndFrameN-1));
            tmmat(i,2) = mean(cellsort_spc(i).tm_trace.subtracted_delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        end
        
        % IEM
        if p.doiem
            % Load
            tmmat(i,5) = mean(cellsort_spc(i).iem_trace.subtracted_delta(p.EndFrame:p.EndFrame+p.EndFrameN-1));
            tmmat(i,4) = mean(cellsort_spc(i).iem_trace.subtracted_delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        end
        
        % Group
        tmmat(i,1) = cellsort_spc(i).group_number;
    end
    
    csvwrite(fullfile(spcpaths.fp_out, spcpaths.tm_csv), tmmat);
end

end