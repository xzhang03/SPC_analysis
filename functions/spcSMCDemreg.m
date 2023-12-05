function spcSMCDemreg(mouse, date, run, varargin)
% spcSMCDemreg applies demonsreg to SMC files

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% File variables
addOptional(p, 'sourcetype', 'smcreg'); % can be smcreg or smc
addOptional(p, 'crop', []); % Crop large FLIM images. Auto detected if source is smcreg
addOptional(p, 'force', false); % Force overwrite or not
addOptional(p, 'savephotons', true); % Output photons
addOptional(p, 'savesmc', true); % Output smcs
addOptional(p, 'outputbinxy', 2); % Output image binxy

% Spatial filter variables
addOptional(p, 'previewlocalnorm', false);
addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric); % Sigma for gaussian fit
addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter

% Registration variables
addOptional(p, 'edges', [0 0 0 0]); % Edge values
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
n = spcpaths.n;

% Reference frames
if ~isempty(p.refrange)
    ind_ref = p.refrange;
else
    ind_ref = 1 : n;
end

%% Check tiffs
% Redo reg SMC or not
if exist(fullfile(spcpaths.fp, sprintf(spcpaths.smcdemreg, 1)), 'file')
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
    if exist(fullfile(spcpaths.fp_out,spcpaths.demregtif_photons), 'file')
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
% Initialize cell
tic;
fprintf('Loading smcs... ');
smccell = cell(n,1);
hwait = waitbar(0, sprintf('Loading smc %i/%i', 1, n));
for i = 1 : n
    if mod(i,5) == 0
        waitbar(i/n, hwait, sprintf('Loading smc %i/%i', i, n));
    end
    loaded = load(fullfile(spcpaths.fp, sprintf(spcpaths.(p.sourcetype), i)), '-mat');
    smccell{i} = loaded.smc;
    if i == 1
        p.crop = loaded.crop;
    end
end
close(hwait);
t = toc;
fprintf('Done. %0.1f s.\n', t);

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
    if iter == 1 && isempty(p.crop) && strcmpi(p.sourcetype, 'smc')
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
    elseif iter == 1 && strcmpi(p.sourcetype, 'smc')
        figure
        imshow(imref,[]);
        rectangle('Position', [p.crop(1) p.crop(2) p.crop(3)-p.crop(1) p.crop(4)-p.crop(2)],...
            'EdgeColor','g', 'LineWidth',2)
    end
    imref = imref(p.crop(2):p.crop(4), p.crop(1):p.crop(3));
    [rows, cols] = size(imref);
    
    %% Apply edge
    if any(p.edges ~= 0)
        % Reference
        imref = imref(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :);

        if mod(size(imref,2), p.binxy) ~= 0 || mod(size(imref,1), p.binxy) ~= 0
            fprintf('Wrong edge values.\n')
            fprintf('Residual in the first two values: %i.\n', mod(size(imref,2),p.binxy));
            fprintf('Residual in the last two values: %i.\n', mod(size(imref,1), p.binxy));
            return;
        end
    end
    
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
    
    %% Register
    tic;
    fprintf('Registering... ');
    hwait = waitbar(0, sprintf('Registering frame%i/%i', 1, n));
    for iframe = 1 : n
        
        % Update wait bar
        if mod(iframe, 5) == 0
            waitbar(iframe/n, hwait, sprintf('Registering frame %i/%i', iframe, n));
        end
        
        % Preprocess frame
        frame = spcSMCphotonframe(smccell{iframe});
        frame_c = frame(p.crop(2):p.crop(4), p.crop(1):p.crop(3));
        frame_ce = frame_c(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :);
        if p.binxy > 1
            frame_ce = binxy(frame_ce, p.binxy);
        end
        frame_cebln = localnormalizecore(frame_ce, p.hp_norm_sigmas, p.medfilt2size);

        % Demons reg
        [D,~] = imregdemons(frame_cebln, imref, [32 16 8 4],...
                    'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitbar',false);
               
        % resize back if necessary
        if p.binxy > 1
            D = imresize(p.binxy * D, p.binxy); % resize to bring back to full size of movie
        end
        D = floor(D);
        
        % Re-embed D-combined into the full-size version if using edges
        if any(p.edges ~= 0)
            D_full = zeros(rows, cols, 2);
            D_full(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :) = D;
            
            % Use the new D_combined for subsequent processing
            D = D_full;
        end
        
        % Reapply crop to make the shift full size
        D_full = zeros(dim(1), dim(2), 2);
        D_full(p.crop(2):p.crop(4), p.crop(1):p.crop(3),:) = D;
        
        % Get D1 and D2
        D1 = D_full(:,:,1);
        D2 = D_full(:,:,2);
        
        % Apply to smc
        smctemp = smccell{iframe};
        for ind = 1 : length(smctemp)
            % Skip empty frames
            if isempty(smctemp(ind).r)
                continue;
            end
            
            % Get r and c
            r = double(smctemp(ind).r);
            c = double(smctemp(ind).c);
            li = r + (c-1) * dim(1);
            
            % Get shifts
            rminus = D2(li);
            cminus = D1(li);
            
            % Apply shifts
            r = r - rminus;
            c = c - cminus;
            
            % Wrap around
            r(r < 1) = r(r < 1) + dim(1);
            r(r > dim(1)) = r(r > dim(1)) - dim(1);
            c(c < 1) = c(c < 1) + dim(2);
            c(c > dim(2)) = c(c > dim(2)) - dim(2);
            
            % Put back
            smctemp(ind).r = uint32(r);
            smctemp(ind).c = uint32(c);
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
hwait = waitbar(0, sprintf('Saving frame %i/%i', 1, n));
for iframe = 1 : n
    % Update wait bar
    if mod(iframe, 5) == 0
        waitbar(iframe/n, hwait, sprintf('Saving frame %i/%i', iframe, n));
    end
        
    % SMC
    if dosmc
        savestruct = struct('smc', smccell{iframe}, 'crop', p.crop); %#ok<NASGU>
        save(fullfile(spcpaths.fp, sprintf(spcpaths.smcdemreg,iframe)), '-struct', 'savestruct');
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
    writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.demregtif_photons), 'double');
end

end

