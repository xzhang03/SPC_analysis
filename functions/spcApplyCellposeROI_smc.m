function spcApplyCellposeROI_smc(mouse, date, run, varargin)
% spcApplySbxROI applies Sbx ROIs to FLIM images

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)
addOptional(p, 'autoc', true); % Automatic cdigit

% Force
addOptional(p, 'force', false); % Force

% Registration variables
addOptional(p, 'sourcetype', 'smcdemreg'); % Input type can be 'smc' or 'smcreg' or 'smcdemreg'

% Binning
addOptional(p, 'binxy', 1);

% Threshold
addOptional(p, 'minarea', 20); % Min number of pixels for a cell (after binning)

% Neuropil
addOptional(p, 'npsize', [16, 4]); % Neuropil size, avoidance size
addOptional(p, 'minnparea', 100); % Minimal number of pixels for neuropil

% Frame info for dff and dt
addOptional(p, 'StartFrame', 1); % Which frame to start
addOptional(p, 'StartFrameN', 12); % How many frames to read for the initial stack
addOptional(p, 'EndFrame', []); % Which frame to begin the end stack. Leave empty to read the last X frames
addOptional(p, 'EndFrameN', 12); % How many frames to read for the end stack

% Time domain: Tm and IEM
addOptional(p, 'tbins', 256); % Time bins
addOptional(p, 'tcycle', 12500); % in ps
addOptional(p, 'T1', 21); % data before this is not considered for tm and iem
addOptional(p, 'T2', 240); % data after this is not considered for tm and iem

% Time domain: IEM99 (using 99th percentile instead of max for IEM peak
% estimates)
addOptional(p, 'iem99', false);
addOptional(p, 'iempercentile', 99);

% Save REF ROI
addOptional(p, 'saverefim', true);

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

% Automate cdigit
if p.autoc
    if ~exist(fullfile(spcpaths.fp, sprintf(spcpaths.sdt_in, 1)), 'file')
        % If can't find the first photon file
        flist = dir(fullfile(spcpaths.fp, sprintf('*.sdt')));
        [cs, ce] = regexp(flist(1).name, 'c\d+.sdt');
        p.cdigit = ce - cs - 4;
        spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
            'slice', p.slice, 'cdigit', p.cdigit);
        fprintf('cdigit updated to %i\n', p.cdigit);
    end
end
    
% Image paths
[smccell, nsections, smccrop] = spcSMCLoad(mouse, date, run, 'spcpaths', spcpaths, 'output', 'cell',...
    'sourcetype',p.sourcetype);

fprintf('Cropping and binning... ')
tic
% Crop
smccell = spcSMCcrop(smccell, smccrop);

% Bin
if p.binxy > 1
    smccell = spcSMCbin(smccell, p.binxy);
end

% Fix endframe
if isempty(p.EndFrame)
    p.EndFrame = nsections - p.EndFrameN + 1;
end

t = toc;
fprintf('Done. %0.1f s.\n', t);

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
trace_ini = struct('raw', zeros(nsections,1), 'neuropil', zeros(nsections,1), ...
    'subtracted', zeros(nsections,1), 'raw_delta', [], 'neuropil_delta', [], 'subtracted_delta', []); % Delta is dff or dt
cellsort_spc = struct('id', [], 'mask', [], 'area', [], 'neuropil', [], 'neuropilarea', [], 'group_number', [],...
    'photon_trace', trace_ini, 'tm_trace', trace_ini, 'tm_start',[], 'iem_trace', trace_ini, 'iem_start',[],...
    'raw_traces', uint16(zeros(p.tbins, nsections)), 'raw_nptraces', uint16(zeros(p.tbins, nsections)),...
    'subtract_traces', uint16(zeros(p.tbins, nsections)));
cellsort_spc = repmat(cellsort_spc, [ncells, 1]);


%% Calculate masks
% Masks
fprintf('Calculating masks and nps... ')
tic

% Neuropil elements
npelbig = strel('disk', p.npsize(1));
npelsmall = strel('disk', p.npsize(2));

% Bin
grinface = binxy(p.grinface, p.binxy) > 0;

% Index for real cells
irealcell = 0;

for i = 1 : ncells
    % Current mask
    maskcurr = masks == i;
    
    % Bin current mask
    maskcurr = binxy(maskcurr, p.binxy) > 0;
    maskarea = sum(maskcurr(:));
    
    if i == 1
        % All ROIs
        allmasks = imdilate(binxy(masks, p.binxy), npelsmall) > 0;

        % Initialize a place to hold all the neuropils
        npall = zeros(size(maskcurr));
    end
    
    % Check area
    if maskarea >= p.minarea
        % Advance real cell counter
        irealcell = irealcell + 1;

        % Neuropil (before binning)
        npcurr = imdilate(maskcurr, npelbig);
        if p.GRIN
            npcurr = npcurr .* double(grinface);
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
                    npcurr = npcurr .* grinface;
                end
                npcurr = npcurr .* ~allmasks;
                nparea = sum(npcurr(:) > 0);
            end
        end

        % Save id
        cellsort_spc(irealcell).id = i;

        % Save mask (a large version)
        cellsort_spc(irealcell).mask = repmat(maskcurr, [1 1 p.tbins]);

        % Save area
        cellsort_spc(irealcell).area = maskarea;

        % Save neuropil (a large version)
        cellsort_spc(irealcell).neuropil = repmat(npcurr > 0, [1 1 p.tbins]);
        npall = npall + npcurr;

        % Save neuropil area
        cellsort_spc(irealcell).neuropilarea = nparea;       
    end
end

% Shrink cellsort
nrealcell = irealcell;
cellsort_spc = cellsort_spc(1:nrealcell);

t = toc;
fprintf('Done. %0.1f s.\n', t);


%% Loop and calculate
% Time resolution
tres = p.tcycle / p.tbins;

% Time vector
ivec = (p.T1 : p.T2)';
tvec = (ivec - ivec(1)) * tres;


% Apply filters
hwait = waitbar(0, 'Extracting traces');
for ii = 1 : nsections
    waitbar(ii/nsections, hwait, sprintf('Extracting traces %i/%i', ii, nsections));
    
    % Get movie
    im = spcDecompress(smccell{ii});
    im = single(im);
    
    for i = 1 : nrealcell
        
        
        % Get trace
        trace = squeeze(sum(sum(im .* cellsort_spc(i).mask,1),2));
        nptrace = squeeze(sum(sum(im .* cellsort_spc(i).neuropil,1),2));
        subtract_trace = trace - nptrace / cellsort_spc(i).neuropilarea * cellsort_spc(i).area;
        subtract_trace(subtract_trace < 0) = 0;
        
        % Load up traces
        cellsort_spc(i).raw_traces(:,ii) = trace;
        cellsort_spc(i).raw_nptraces(:,ii) = nptrace;
        cellsort_spc(i).subtract_traces(:,ii) = subtract_trace;
        
        % Trim
        trace = trace(ivec);
        nptrace = nptrace(ivec);
        subtract_trace = subtract_trace(ivec);
        
        % Photon trace value
        cellsort_spc(i).photon_trace.raw(ii) = sum(trace);
        cellsort_spc(i).photon_trace.neuropil(ii) = sum(nptrace);
        cellsort_spc(i).photon_trace.subtracted(ii) = cellsort_spc(i).photon_trace.raw(ii)...
            - cellsort_spc(i).photon_trace.neuropil(ii);
        
        % TM trace value
        cellsort_spc(i).tm_trace.raw(ii) = sum(trace .* tvec) / cellsort_spc(i).photon_trace.raw(ii);
        cellsort_spc(i).tm_trace.neuropil(ii) = sum(nptrace .* tvec) / cellsort_spc(i).photon_trace.neuropil(ii);
        cellsort_spc(i).tm_trace.subtracted(ii) = sum(subtract_trace .* tvec) / sum(subtract_trace);
        
        % IEM trace value
        if p.iem99
            cellsort_spc(i).iem_trace.raw(ii) = cellsort_spc(i).photon_trace.raw(ii) / prctile(trace, p.iempercentile) * tres;
            cellsort_spc(i).iem_trace.neuropil(ii) = cellsort_spc(i).photon_trace.neuropil(ii) / prctile(trace, p.iempercentile) * tres;
            cellsort_spc(i).iem_trace.subtracted(ii) = sum(subtract_trace) / prctile(trace, p.iempercentile) * tres;
        else
            cellsort_spc(i).iem_trace.raw(ii) = cellsort_spc(i).photon_trace.raw(ii) / max(trace) * tres;
            cellsort_spc(i).iem_trace.neuropil(ii) = cellsort_spc(i).photon_trace.neuropil(ii) / max(nptrace) * tres;
            cellsort_spc(i).iem_trace.subtracted(ii) = sum(subtract_trace) / max(subtract_trace) * tres;
        end
    end
end
close(hwait)

%% Calculate changes
for i = 1 : nrealcell
    %% Photons
    % Photon traces
    dff_trace = cellsort_spc(i).photon_trace.raw ...
        / mean(cellsort_spc(i).photon_trace.raw(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
    cellsort_spc(i).photon_trace.raw_delta = dff_trace;

    % dff photon np trace
    dff_np_trace = cellsort_spc(i).tm_trace.neuropil ...
        / mean(cellsort_spc(i).tm_trace.neuropil(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
    cellsort_spc(i).photon_trace.neuropil_delta = dff_np_trace;

    % dff photon subtracted trace
    dff_subtracted_trace = cellsort_spc(i).photon_trace.subtracted ...
        / mean(cellsort_spc(i).photon_trace.subtracted(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
    cellsort_spc(i).photon_trace.subtracted_delta = dff_subtracted_trace;
    
    %% TM
    % dt tm trace
    dt_trace = cellsort_spc(i).tm_trace.raw...
        - mean(cellsort_spc(i).tm_trace.raw(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    cellsort_spc(i).tm_trace.raw_delta = dt_trace;

    % dt tm np trace
    dt_np_trace = cellsort_spc(i).tm_trace.neuropil...
        - mean(cellsort_spc(i).tm_trace.neuropil(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    cellsort_spc(i).tm_trace.neuropil_delta = dt_np_trace;

    % dt tm subtracted trace
    dt_subtracted_trace = cellsort_spc(i).tm_trace.subtracted...
        - mean(cellsort_spc(i).tm_trace.subtracted(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    cellsort_spc(i).tm_trace.subtracted_delta = dt_subtracted_trace;

    % Tm start
    cellsort_spc(i).tm_start = mean(cellsort_spc(i).tm_trace.raw(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    
    %% IEM
    % dt iem trace
    dt_trace = cellsort_spc(i).iem_trace.raw...
        - mean(cellsort_spc(i).iem_trace.raw(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    cellsort_spc(i).iem_trace.raw_delta = dt_trace;

    % dt iem np trace
    dt_np_trace = cellsort_spc(i).iem_trace.neuropil...
        - mean(cellsort_spc(i).iem_trace.neuropil(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    cellsort_spc(i).iem_trace.neuropil_delta = dt_np_trace;

    % dt iem subtracted trace
    dt_subtracted_trace = cellsort_spc(i).iem_trace.subtracted...
        - mean(cellsort_spc(i).iem_trace.subtracted(p.StartFrame : p.StartFrame+p.StartFrameN-1));
    cellsort_spc(i).iem_trace.subtracted_delta = dt_subtracted_trace;

    % IEM start
    cellsort_spc(i).iem_start = mean(cellsort_spc(i).iem_trace.raw(p.StartFrame : p.StartFrame+p.StartFrameN-1));
end

%% Save to file
% Expand npall
npall = imresize(npall, p.binxy);
if size(npall,1) < size(masks,1)
    npall = cat(1, npall, zeros(size(masks,1)-size(npall,1), size(npall,2)));
end
if size(npall,2) < size(masks,2)
    npall = cat(2, npall, zeros(size(npall,1), size(masks,2)-size(npall,2)));
end

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
    ROIs_bin = npall > 0.5;
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

end