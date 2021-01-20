function spcApplySbxROI(mouse, date, run, varargin)
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

% Sbx variables
addOptional(p, 'sbxrun', []); % Sbx run number
addOptional(p, 'checkwarp', false); % Check warp

% Cropping
addOptional(p, 'margin', 0);

% Binning
addOptional(p, 'binxy', 1);

% Pick datasets
addOptional(p, 'dophotons', true); % Photon count
addOptional(p, 'dotm', true); % Mean lifetime
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

% Get sbx path
sbxpath = sbxPath(mouse, date, p.sbxrun, 'signals', 'server', p.server);

% Image paths
% Switch types
switch p.sourcetype
    case 'raw'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.tif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.tif_tm);
    case 'registered'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.regtif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.regtif_tm);
    case 'warped'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.warptif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.warptif_tm);
    case 'demonsreg'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.demregtif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.demregtif_tm);
end

%% Check warp
if p.checkwarp
    % Load sbx sample
    % Folders and files
    datemouse = sprintf('%s_%s', date, mouse);
    datemouserun = sprintf('%s_%s_run%i', date, mouse, p.sbxrun);
    sbxrefpath = sprintf('\\\\%s\\data\\2p\\%s\\%s\\%s\\%s\\',...
    p.server, p.user, mouse, datemouse, datemouserun);
    sbxrefname = sprintf('%s_%s_run%i_ref.tif', date, mouse, p.sbxrun);
    
    % Load
    sbxsample = readtiff(fullfile(sbxrefpath, sbxrefname));
    
    % Load spc photon data
    im_photon = readtiff(fp_photon);
    
    % Make spc sample
    spcsample = mean(im_photon, 3);
    
    % Crop
    [sbxCropped, spcCropped] = xregCropping(sbxsample, spcsample, p.margin);
    
    figure
    imshowpair(sbxCropped, spcCropped);
    title('Check warping')
end

%% Load data
% Load the corresponding data
sigstruct = load(sbxpath, '-mat');

% Basic paramteres
ncells = length(sigstruct.cellsort) - 1;


%% Initialize
% datastruct
if isfield(sigstruct, 'cellsort_spc') && ~p.force
    cellsort_spc = sigstruct.cellsort_spc;
    
    % Redo?
    if ~isempty(cellsort_spc(end).photon_trace.raw) && p.dophotons
        p.dophotons = questdlg('Photon trace already exists. Redo?', 'Redo photon trace', 'Yes', 'No', 'No');
        p.dophotons = strcmp(p.dophotons, 'Yes');
    end
    
    if ~isempty(cellsort_spc(end).tm_trace.raw) && p.dotm
        p.dotm = questdlg('Tm trace already exists. Redo?', 'Redo Tm trace', 'Yes', 'No', 'No');
        p.dotm = strcmp(p.dotm, 'Yes');
    end
else
    trace_ini = struct('raw', [], 'delta', []); % Delta is dff or dt
    cellsort_spc = struct('mask', [], 'area', [], 'group_number', [],...
        'photon_trace', trace_ini, 'tm_trace', trace_ini);
    cellsort_spc = repmat(cellsort_spc, [ncells, 1]);
end

%% Photons and Tm
% Load photon data
if p.dophotons && ~exist('im_photon', 'var')
    % Load warped
    [im_photon, ~] = readtiff(fp_photon);
end

% Load tm data
if p.dotm
    % Load warped
    [im_tm, ~] = readtiff(fp_tm);
end

% Get the number of sections
if exist('im_photon', 'var')
    nsections = size(im_photon, 3);
elseif exist('im_tm', 'var')
    nsections = size(im_tm, 3);
else
    nsections = length(cellsort_spc(1).tm_trace.raw);
end

% Fix endframe
if isempty(p.EndFrame)
    p.EndFrame = nsections - p.EndFrameN + 1;
end

% Get sbx and spc canvas samples for cropping
sbxcanvas = sigstruct.movROI;
if p.dophotons
    spccanvas = mean(im_photon, 3);
elseif p.dotm
    spccanvas = mean(im_tm, 3);
end

% Crop photon data
if p.dophotons
    % First frame
    [~, im_photon_crop] = xregCropping(sbxcanvas, im_photon(:,:,1), p.margin);
    
    % Initialize
    im_photon_crop = repmat(im_photon_crop, [1 1 nsections]);
    
    % Loop and crop
    for i = 2 : nsections
        [~, im_photon_crop(:,:,i)] = xregCropping(sbxcanvas, im_photon(:,:,i), p.margin);
    end
    
end

% Crop tm data
if p.dotm
    % First frame
    [~, im_tm_crop] = xregCropping(sbxcanvas, im_tm(:,:,1), p.margin);
    
    % Initialize
    im_tm_crop = repmat(im_tm_crop, [1 1 nsections]);
    
    % Loop and crop
    for i = 2 : nsections
        [~, im_tm_crop(:,:,i)] = xregCropping(sbxcanvas, im_tm(:,:,i), p.margin);
    end
end

% Bin
if p.dophotons && p.binxy > 1
    im_photon_crop = binxy(im_photon_crop, p.binxy);
end
if p.dotm && p.binxy > 1
    im_tm_crop = binxy(im_tm_crop, p.binxy);
end

if p.dotm || p.dophotons
    % Apply filters
    hwait = waitbar(0, 'Applying masks');
    for i = 1 : ncells
        if mod(i, 10) == 0
            waitbar(i/ncells);
        end

        % Crop filter
        [maskcurr, ~] = xregCropping(sigstruct.cellsort(i).mask, spccanvas, p.margin);

        % Bin filter if needed
        if p.binxy > 1
            maskcurr = binxy(maskcurr, p.binxy);
        end

        % Save mask
        cellsort_spc(i).mask = maskcurr;

        % Save group
        cellsort_spc(i).group_number = sigstruct.cellsort(i).group_number;

        % Area
        currarea = sum(maskcurr(:));

        % Save area
        cellsort_spc(i).area = currarea;

        % Replicate filter for sections
        mask_stack = repmat(maskcurr, [1 1 nsections]);

        % Get photon counts
        if p.dophotons
            % raw photon trace
            photon_trace =...
                squeeze(nansum(nansum(im_photon_crop .* mask_stack, 1), 2)) / currarea;
            cellsort_spc(i).photon_trace.raw = photon_trace;

            % dff photon trace
            dff_trace =...
                photon_trace / mean(photon_trace(p.StartFrame : p.StartFrame+p.StartFrameN)) - 1;
            cellsort_spc(i).photon_trace.delta = dff_trace;
        end
        if p.dotm
            % raw tm trace
            tm_trace =...
                squeeze(nansum(nansum(im_tm_crop .* mask_stack, 1), 2)) / currarea;
            cellsort_spc(i).tm_trace.raw = tm_trace;

            % dt tm trace
            dt_trace =...
                tm_trace - mean(tm_trace(p.StartFrame : p.StartFrame+p.StartFrameN));
            cellsort_spc(i).tm_trace.delta = dt_trace;
        end
    end
    close(hwait)
end
%% Plot
if p.plotscatter
    % Initialize
    photonsvec = nan(ncells,1);
    tmvec = nan(ncells,1);
    groupvec = nan(ncells,1);
    
    for i = 1 : ncells
        % Load
        photonsvec(i) =...
            mean(cellsort_spc(i).photon_trace.delta(p.EndFrame:p.EndFrame+p.EndFrameN-1))...
            - mean(cellsort_spc(i).photon_trace.delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        
        tmvec(i) =...
            mean(cellsort_spc(i).tm_trace.delta(p.EndFrame:p.EndFrame+p.EndFrameN-1))...
            - mean(cellsort_spc(i).tm_trace.delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        
        groupvec(i) = cellsort_spc(i).group_number;
    end

    figure
    gscatter(photonsvec, tmvec, groupvec);
end

%% Save to file
% Save structure
sigstruct.cellsort_spc = cellsort_spc;

% Save file
save(sbxpath, '-struct', 'sigstruct', '-v7.3');

if p.saverefim
    % Initialize
    ROIframe = zeros(size(im_photon_crop,1),size(im_photon_crop,2));
    
    for i = 1 : ncells
        ROIframe = ROIframe + cellsort_spc(i).mask;
    end
        
    % Make rgb
    rgb = repmat(mat2gray(mean(im_photon_crop, 3)),[1 1 3]);
    rgb(:,:,1) = (ROIframe > 0) * 0.5;
    rgb(:,:,3) = 0;
    rgb = imresize(rgb, p.binxy);
    
    % show
    figure
    imshow(rgb)

    % Save figure
    saveas(gcf, fullfile(spcpaths.fp_out, spcpaths.run_ROI_ref));
end

if p.savetmcsv
    % Iinitialize
    tmmat = nan(ncells,3);
    
    for i = 1 : ncells
        % Load
        tmmat(i,3) = mean(cellsort_spc(i).tm_trace.delta(p.EndFrame:p.EndFrame+p.EndFrameN-1));
        tmmat(i,2) = mean(cellsort_spc(i).tm_trace.delta(p.StartFrame:p.StartFrame+p.StartFrameN-1));
        
        % Group
        tmmat(i,1) = cellsort_spc(i).group_number;
    end
    
    csvwrite(fullfile(spcpaths.fp_out, spcpaths.tm_csv), tmmat);
end

end