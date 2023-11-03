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

% Time domain: tm and iem
addOptional(p, 'tbins', 256); % Time bins
addOptional(p, 'tcycle', 12500); % in ps
addOptional(p, 'T1', 21); % data before this is not considered for tm and iem
addOptional(p, 'T2', 240); % data after this is not considered for tm and iem
addOptional(p, 'threshold', 5); % Peak has to be above this number for tm and iem considerations
addOptional(p, 'bin_tm', 5); % Square binning for tm calculation. Edge 2x+1
addOptional(p, 'bin_iem', 5); % Square binning for tm calculation. Edge 2x+1

% Resize (Cropping at loading stage)
addOptional(p, 'resizedim', [512 1250]);

% Compress
addOptional(p, 'compress', true); % First time it will make compression file, later reading it

% IRF
addOptional(p, 'deconvforiem', false); % If set true it takes 10x as long
addOptional(p, 'irf', [882; 6176; 15000; 15000; 3529]);

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
hwait = waitbar(0);
for ind = spcpaths.cinds
    % Load movie
    tic
    if p.compress
        usesmc = exist(fullfile(spcpaths.fp, sprintf(spcpaths.smc,ind)), 'file');
    else
        usesmc = false;
    end
    
    if usesmc
        % Load sparse matrix compression file
        smc = load(fullfile(spcpaths.fp, sprintf(spcpaths.smc,ind)), "-mat");
        smc = smc.smc;
        mov = spcDecompress(smc);
    else
        mov = spcLoadsdt(fullfile(spcpaths.fp, sprintf(spcpaths.sdt_in,ind)), 'crop', p.resizedim);

        if p.compress
            % Save smc for future
            smc = spcCompress(mov);
            save(fullfile(spcpaths.fp, sprintf(spcpaths.smc,ind)), 'smc', '-v7.3');
        end
    end
    t = toc;
    fprintf('Loaded Frame %i in %0.2f s.', ind, t);
    
    % Crop
    if ind == 1 && isempty(p.crop)
        figure
        imshow(sum(mov,3),[]);
        p.crop = round(wait(imrect()));
        p.crop(3) = min(p.crop(1) + p.crop(3), size(mov,2));
        p.crop(1) = max(p.crop(1),1);
        p.crop(4) = min(p.crop(2) + p.crop(4), size(mov,1));
        p.crop(2) = max(p.crop(2),1);
        close(gcf)
    end
    mov = mov(p.crop(2):p.crop(4), p.crop(1):p.crop(3), :);
    sz = size(mov);
    
    % Initialize
    if dophotons
        curr_photon = sum(mov,3);
    end
    if dotm
        curr_tm = zeros(sz(1), sz(2));
    end
    if doiem
        curr_iem = zeros(sz(1), sz(2));
    end
    
    
    % Time processing
    for r = 1 : sz(1)
        if mod(r,100) == 0
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
                    curr_tm(r,c) = nan;
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
                    curr_iem(r,c) = sum(trace) / max(trace) * tres;
                else
                    curr_iem(r,c) = nan;
                end
            end
        end
    end
    t = toc;
    fprintf(' Processing %0.2f s.\n', t);
    
    % Repmat
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
end
close(hwait)

%% Save
if dophotons
    % Write
    writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.tif_photons), 'double');
end
if dotm
    % Write
    writetiff(im_tm, fullfile(spcpaths.fp_out,spcpaths.tif_tm), 'double');
end
if doiem
    % Write
    writetiff(im_iem, fullfile(spcpaths.fp_out,spcpaths.tif_iem), 'double');
end

% Parameters
save(fullfile(spcpaths.fp_out, spcpaths.loadparams), 'p', '-v7.3');

end