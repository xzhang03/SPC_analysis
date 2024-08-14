function spcTiff_sdt(mouse, date, run, varargin)
% spcTiff converts raw sdt files to tiff files

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice

addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)
addOptional(p, 'autoc', false); % Automate cdigit - experimental

% File variables
addOptional(p, 'photons', true); % Make photons file or not
addOptional(p, 'tm', true); % Make tm file or not
addOptional(p, 'iem', true); % Make iem files or not
addOptional(p, 'force', false); % Force overwrite or not

% Crop
addOptional(p, 'crop', []);

% Bin xy
addOptional(p, 'binxy', 1); % This bin determins the spatial resolution of the frame. 
                            % By default, the SPCimage analysis pipeline
                            % always have a binxy of 1. But that's slow
                            % for IEM so now you can bin the frame before
                            % the tau bin. Note that a binxy of 2 an a
                            % bin_tm of 3, means an overall bin of 6.


% Time domain: tm and iem
addOptional(p, 'tbins', 256); % Time bins
addOptional(p, 'tcycle', 12500); % in ps
addOptional(p, 'T1', 21); % data before this is not considered for tm and iem
addOptional(p, 'T2', 240); % data after this is not considered for tm and iem
addOptional(p, 'threshold', 5); % Peak has to be above this number for tm and iem considerations
addOptional(p, 'bin_tm', 5); % Square binning for tm calculation. Edge 2x+1. This is spatial binning on top of binxy.
addOptional(p, 'bin_iem', 5); % Square binning for tm calculation. Edge 2x+1. This is spatial on top of binxy.

% Time domain: IEM99 (using 99th percentile instead of max for IEM peak
% estimates)
addOptional(p, 'iem99', false);
addOptional(p, 'iempercentile', 99);

% Resize (Cropping at loading stage)
addOptional(p, 'resizedim', [512 1250]);

% Compress
addOptional(p, 'compress', true); % First time it will make compression file, later reading it
addOptional(p, 'smctype', ''); % Leave empty to use sdt or smc. Specify for smcreg or smcdemreg

% IRF
addOptional(p, 'deconvforiem', false); % If set true it takes 10x as long
addOptional(p, 'irf', [882; 6176; 15000; 15000; 3529]);

% Live preview
addOptional(p, 'livepreview', true);
addOptional(p, 'previewim', 'photons');

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
        case 'MP'
            p.user = 'marta';
    end
end

%% Time vectors
% Time resolution
tres = p.tcycle / p.tbins;

% Time vector
ivec = (p.T1 : p.T2)';
tvec = (ivec - p.T1) * tres;

%% IO
% Get paths
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

%% Check if redo
% Do photon data?
if p.photons
    if exist(fullfile(spcpaths.fp_out,spcpaths.tif_photons), 'file')
        % File already exists
        if p.force
            % Force redo
            dophotons = true;
        elseif input('Photon tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
            % Ask redo
            dophotons = true;
        else
            % No redo
            dophotons = false;
        end
    else
        % No file
        dophotons = true;
    end
else
    % Do nothing
    dophotons = false;
end

% Do Tm data?
if p.tm
    if exist(fullfile(spcpaths.fp_out,spcpaths.tif_tm), 'file')
        % File already exists
        if p.force
            % Force redo
            dotm = true;
        elseif input('Tm tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
            % Ask redo
            dotm = true;
        else
            % No redo
            dotm = false;
        end
    else
        % No file
        dotm = true;
    end
else
    % Do nothing
    dotm = false;
end

% Do IEM data?
if p.iem
    if exist(fullfile(spcpaths.fp_out,spcpaths.tif_iem), 'file')
        % File already exists
        if p.force
            % Force redo
            doiem = true;
        elseif input('IEM tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
            % Ask redo
            doiem = true;
        else
            % No redo
            doiem = false;
        end
    else
        % No file
        doiem = true;
    end
else
    % Do nothing
    doiem = false;
end

if ~dophotons && ~dotm && ~doiem
    disp('Nothing to do.')
    return;
end


%% Load
hwait = waitbar(0, sprintf('Processing %s %s run%i: %i/%i r%i', mouse, date, run, 1, spcpaths.n, 0));
for ind = spcpaths.cinds
    %% Load movie
    tic
    if isempty(p.smctype)
        if p.compress
            usesmc = exist(fullfile(spcpaths.fp, sprintf(spcpaths.smc,ind)), 'file');
            p.smctype = 'smc';
        else
            usesmc = false;
        end
    else
        usesmc = true;
    end
    
    if usesmc
        % Load sparse matrix compression file
        loaded = load(fullfile(spcpaths.fp, sprintf(spcpaths.(p.smctype),ind)), "-mat");
        if ind == 1 && isfield(loaded, 'crop')
            if ~isempty(loaded.crop)
                p.crop = loaded.crop;
            end
        end
        
        smc = loaded.smc;
        mov = spcDecompress(smc);
    else
        % Load sdt
        mov = spcLoadsdt(fullfile(spcpaths.fp, sprintf(spcpaths.sdt_in,ind)), 'crop', p.resizedim);

        if p.compress
            % Save smc for future
            smc = spcCompress(mov);
            save(fullfile(spcpaths.fp, sprintf(spcpaths.smc,ind)), 'smc', '-v7.3');
        end
    end
    t = toc;
    fprintf('Loaded Frame %i in %0.2f s.', ind, t);
    tic;
    
    %% Crop
    if ind == 1 && isempty(p.crop)
        figure
        imshow(sum(mov,3),[]);
        p.crop = round(wait(imrect()));
        p.crop(3) = min(p.crop(1) + p.crop(3), size(mov,2));
        p.crop(1) = max(p.crop(1),1);
        p.crop(4) = min(p.crop(2) + p.crop(4), size(mov,1));
        p.crop(2) = max(p.crop(2),1);
        disp(p.crop);
        close(gcf)
    elseif ind == 1
        figure
        imshow(sum(mov,3),[]);
        rectangle('Position', [p.crop(1) p.crop(2) p.crop(3)-p.crop(1) p.crop(4)-p.crop(2)],...
            'EdgeColor','g', 'LineWidth',2)
    end
    mov = mov(p.crop(2):p.crop(4), p.crop(1):p.crop(3), :);
    
    %% XY bin
    if p.binxy > 1
        mov = binxy(mov, p.binxy) * (p.binxy^2);
    end
    sz = size(mov);
    
    %% Initialize
    if dophotons
        curr_photon = sum(mov,3);
    end
    if dotm
        curr_tm = zeros(sz(1), sz(2));
    end
    if doiem
        curr_iem = zeros(sz(1), sz(2));
    end
    
    
    %% Time processing
    for r = 1 : sz(1)
        if mod(r,50) == 0
            waitbar(ind/spcpaths.n, hwait, sprintf('Processing %s %s run%i: %i/%i r%i', mouse, date, run, ind, spcpaths.n, r));
        end
        
        for c = 1 : sz(2)
            % Tm calculation
            if dotm
                trace = spcMovtrace(mov, r, c, p.bin_tm);
                trace = trace(ivec);
                if any(trace > p.threshold)
                    curr_tm(r,c) = sum(trace .* tvec) / sum(trace);
                else
                    curr_tm(r,c) = 0;
                end  
            end
            
            % IEM calculation
            if doiem
                % Resample traces if needed
                if p.bin_iem ~= p.bin_tm
                    trace = spcMovtrace(mov, r, c, p.bin_iem);
                    trace = trace(ivec);
                end
                if any(trace > p.threshold)
                    if p.deconvforiem
                        trace = deconvlucy(trace,p.irf);
                    end
                    
                    if p.iem99
                        % IEM 99
                        curr_iem(r,c) = sum(trace) / prctile(trace, p.iempercentile) * tres;
                    else
                        curr_iem(r,c) = sum(trace) / max(trace) * tres;
                    end
                else
                    curr_iem(r,c) = 0;
                end
            end
        end
    end
    t = toc;
    fprintf(' Processing %0.2f s.\n', t);
    
    %% Repmat and save
    if ind == 1
        % Initialize if the first frame
        if dophotons
            im_photon = repmat(curr_photon, [1 1 length(spcpaths.cinds)]);
        end
        if dotm
            im_tm = repmat(curr_tm, [1 1 length(spcpaths.cinds)]);
        end
        if doiem
            im_iem = repmat(curr_iem, [1 1 length(spcpaths.cinds)]);
        end
    else
        if dophotons
            im_photon(:,:,ind) = curr_photon;
        end
        if dotm
            im_tm(:,:,ind) = curr_tm;
        end
        if doiem
            im_iem(:,:,ind) = curr_iem;
        end
    end
    
    %% Preview
    if p.livepreview
        if ind == 1
            % New figure
            figure;
            
            % Show
            switch p.previewim
                case 'photons'
                    him = imshow(curr_photon, []);
                case 'tm'
                    him = imshow(curr_tm, []);
                case 'iem'
                    him = imshow(curr_iem, []);
            end
        else
            % Update
            switch p.previewim
                case 'photons'
                    him.CData = curr_photon;
                    him.Parent.CLim(2) = max(curr_photon(:));
                case 'tm'
                    him.CData = curr_tm;
                    him.Parent.CLim(2) = max(curr_tm(:));
                case 'iem'
                    him.CData = curr_iem;
                    him.Parent.CLim(2) = max(curr_iem(:));
            end
        end
        
        title(sprintf('Preview %i/%i', ind, spcpaths.n))
    end
end
close(hwait)

%% Save
if dophotons
    % Write
    switch p.smctype
        case 'smcreg'
            writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.regtif_photons), 'double');
        case 'smcdemreg'
            writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.demregtif_photons), 'double');
        otherwise
            writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.tif_photons), 'double');
    end
end
if dotm
    % Write
    switch p.smctype
        case 'smcreg'
            writetiff(im_tm, fullfile(spcpaths.fp_out,spcpaths.regtif_tm), 'double');
        case 'smcdemreg'
            writetiff(im_tm, fullfile(spcpaths.fp_out,spcpaths.demregtif_tm), 'double');
        otherwise
            writetiff(im_tm, fullfile(spcpaths.fp_out,spcpaths.tif_tm), 'double');
    end
end
if doiem
    % Write
    if p.iem99
        % IEM99
        switch p.smctype
            case 'smcreg'
                writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.regtif_iem99), 'double');
            case 'smcdemreg'
                writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.demregtif_iem99), 'double');
            otherwise
                writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.tif_iem99), 'double');
        end
    else
        % Regular iem
        switch p.smctype
            case 'smcreg'
                writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.regtif_iem), 'double');
            case 'smcdemreg'
                writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.demregtif_iem), 'double');
            otherwise
                writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.tif_iem), 'double');
        end
    end
end

% Parameters
save(fullfile(spcpaths.fp_out, spcpaths.loadparams), 'p', '-v7.3');

end