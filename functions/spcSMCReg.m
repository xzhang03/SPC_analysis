function spcSMCReg(mouse, date, run, varargin)
% spcSMCReg register from SMC compressed data. This is XY reg only

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% File variables
addOptional(p, 'crop', []); % Crop
addOptional(p, 'force', false); % Force overwrite or not
addOptional(p, 'savephotons', true); % Output photons
addOptional(p, 'outputbinxy', 2); % Output image binxy

% Spatial filter variables
addOptional(p, 'previewlocalnorm', true);
addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric); % Sigma for gaussian fit
addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter

% Registration variables
addOptional(p, 'binxy', 1); % Binning
addOptional(p, 'refrange', []); % Specify first and last frame numbers, leave empty to us all
addOptional(p, 'iterations', 1); % Repeated iterations

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

%% IO
% Get paths
spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);

%% Check tiffs
% Redo reg SMC or not
if exist(fullfile(spcpaths.fp, sprintf(spcpaths.smcreg, 1)), 'file')
    % File already exists
    if p.force
        % Force redo
        dosmc = true;
    elseif input('SMC reg files already exists, redo? (1 = yes, 0 = no): ') == 1
        % Ask redo
        dosmc = true;
    else
        % No redo
        dosmc = false;
    end
else
    % No file
    dosmc = true;
end
    
% Redo photons or not
if p.savephotons
    if exist(fullfile(spcpaths.fp_out,spcpaths.regtif_photons), 'file')
        % File already exists
        if p.force
            % Force redo
            dophotons = true;
        elseif input('Photon reg tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
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
    dophotons = false;
end

if ~dosmc && ~dophotons
    disp('Nothing to do');
    return;
end

%% Reading
[smccell, n, smccrop] = spcSMCLoad(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit, 'sourcetype', 'smc', 'output', 'cell');
if ~isempty(smccrop)
    p.crop = smccrop;
end

% Reference frames
if ~isempty(p.refrange)
    ind_ref = p.refrange;
else
    ind_ref = 1 : n;
end

%% Loop through iterations
for iter = 1 : p.iterations
    
    %% Constructing reference
    tic;
    fprintf('Making ref... ');
    fref = true; % First ref frame flag
    for iref = ind_ref
        imref_t = spcSMCphotonframe(smccell{iref});
        if fref
            % First frame
            imref = imref_t;
            fref = false;
        else
            imref = imref + imref_t;
        end
    end
    t = toc;
    fprintf('Done. %0.1f s.\n', t);
    
    %% Crop
    dim = size(imref);
    if iter == 1 && isempty(p.crop)
        figure
        imshow(sum(imref,3),[]);
        title('Cropping.')
        p.crop = round(wait(imrect()));
        p.crop(3) = min(p.crop(1) + p.crop(3), dim(2));
        p.crop(1) = max(p.crop(1),1);
        p.crop(4) = min(p.crop(2) + p.crop(4), dim(1));
        p.crop(2) = max(p.crop(2),1);
        disp(p.crop);
        close(gcf)
    elseif iter == 1
        figure
        imshow(imref,[]);
        rectangle('Position', [p.crop(1) p.crop(2) p.crop(3)-p.crop(1) p.crop(4)-p.crop(2)],...
            'EdgeColor','g', 'LineWidth',2)
    end
    imref = imref(p.crop(2):p.crop(4), p.crop(1):p.crop(3));
    
    %% Bin ref
    if p.binxy > 1
        imref = binxy(imref, p.binxy);
    end
    
    %% Local norm preview
    if iter == 1 && p.previewlocalnorm
        lnparas = localnormalize_ui('im', imref, 'gausssizes', p.hp_norm_sigmas, 'parametermode', true);
        p.hp_norm_sigmas = lnparas([2,1]);
    end
    
    %% Local normalize ref
    imref = localnormalizecore(imref, p.hp_norm_sigmas, p.medfilt2size);
    figure
    imshow(imref,[]);
    title(sprintf('Iteration %i', iter));
    
    %% Get a focus area
    if iter == 1
        % Get a focus area for calculating shifts
        figure
        imshow(imref,[]);
        title('Focus area.')
        h = imrect;
        regfocus = wait(h);
        regfocus = round(regfocus);
        close(gcf)
        regfocus2 = [regfocus(2), regfocus(2) + regfocus(4) - 1, regfocus(1), regfocus(1) + regfocus(3) - 1];
    end
    imref = imref(regfocus2(1) : regfocus2(2), regfocus2(3) : regfocus2(4));
    
    %% Register
    tic;
    fprintf('Registering... ');
    hwait = waitbar(0, sprintf('Registering frame %i/%i', 1, n));
    for iframe = 1 : n
        % Update wait bar
        if mod(iframe, 5) == 0
            waitbar(iframe/n, hwait, sprintf('Registering frame %i/%i', iframe, n));
        end
        
        % Preprocess frame
        frame = spcSMCphotonframe(smccell{iframe});
        frame_c = frame(p.crop(2):p.crop(4), p.crop(1):p.crop(3));
        if p.binxy > 1
            frame_c = binxy(frame_c, p.binxy);
        end
        frame_cln = localnormalizecore(frame_c, p.hp_norm_sigmas, p.medfilt2size);
        frame_clnf = frame_cln(regfocus2(1) : regfocus2(2), regfocus2(3) : regfocus2(4));
        
        % Register
        [xy_shifts,~] = stackRegisterMA_RR(frame_clnf, imref, 100, [], 0);
        
        % Apply shifts to smc
        rplus = round(xy_shifts(3) * p.binxy);
        cplus = round(xy_shifts(4) * p.binxy);
        
        % smc
        smctemp = smccell{iframe};
        for ind = 1 : length(smctemp)
            % Skip empty frames
            if isempty(smctemp(ind).r)
                continue;
            end
            
            % Shift
            r = uint16(smctemp(ind).r + rplus);
            c = uint16(smctemp(ind).c + cplus);
            
            % Wrap around
            r(r < 1) = r(r < 1) + dim(1);
            r(r > dim(1)) = r(r > dim(1)) - dim(1);
            c(c < 1) = c(c < 1) + dim(2);
            c(c > dim(2)) = c(c > dim(2)) - dim(2);
            
            % Put back
            smctemp(ind).r = r;
            smctemp(ind).c = c;
        end
        smccell{iframe} = smctemp;
    end
    close(hwait);
    t = toc;
    fprintf('Done. %0.1f s.\n', t);
end

%% Save files
tic;
fprintf('Saving smcs... ');
hwait = waitbar(0, sprintf('Saving frame%i/%i', 1, n));
for iframe = 1 : n
    % Update wait bar
    if mod(iframe, 5) == 0
        waitbar(iframe/n, hwait, sprintf('Saving frame %i/%i', iframe, n));
    end
        
    % SMC
    if dosmc
        savestruct = struct('smc', smccell{iframe}, 'crop', p.crop); %#ok<NASGU>
        save(fullfile(spcpaths.fp, sprintf(spcpaths.smcreg,iframe)), '-struct', 'savestruct');
    end
    
    % Decompress
    if dophotons
        % Decompress
        curr_photon = spcSMCphotonframe(smccell{iframe});
        
        % Crop
        curr_photon = curr_photon(p.crop(2):p.crop(4), p.crop(1):p.crop(3));
        
        % XY bin
        if p.outputbinxy > 1
            curr_photon = binxy(curr_photon, p.outputbinxy);
        end
        
        if iframe == 1
            im_photon = repmat(curr_photon, [1 1 n]);
        else
            im_photon(:,:,iframe) = curr_photon;
        end
    end
    
end
close(hwait);
t = toc;
fprintf('Done. %0.1f s.\n', t);

if dophotons
    % Write
    writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.regtif_photons), 'double');
end
end