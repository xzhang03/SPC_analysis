function spcCellposePrep(mouse, date, run, varargin)
% spcTiffCellposePrep preps the dataset for cellpose segmentation

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% File variables
addOptional(p, 'sourcetype', 'registered'); % Can be 'raw', 'reigstered', 'warped'
addOptional(p, 'force', false); % Force overwrite or not

% Spatial filter variables
addOptional(p, 'usemean', true); % Use mean instead of median (better for low photon counts)
addOptional(p, 'uselocalnorm', true);
addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric); % Sigma for gaussian fit
addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter

% Show
addOptional(p, 'saveraw', true); % Save a copy without the filtering in case that's better
addOptional(p, 'showresult', true);

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
% Redo photons or not
if exist(fullfile(spcpaths.fp_out,spcpaths.cp_toseg), 'file')
    % File already exists
    if p.force
        % Force redo
        dophotons = true;
    elseif input('Cellpose-to-seg tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
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

%% Load
if dophotons
    % Switch types
    switch p.sourcetype
        case 'raw'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.tif_photons);
            im_photon = readtiff(fp_photon);
        case 'registered'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.regtif_photons);
            im_photon = readtiff(fp_photon);
        case 'warped'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.warptif_photons);
            im_photon = readtiff(fp_photon);
        case 'demonsreg'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.demregtif_photons);
            im_photon = readtiff(fp_photon);
        otherwise
            [smcstruct, n, smccrop] = spcSMCLoad(mouse, date, run, 'server', p.server, 'user', p.user,...
                'slice', p.slice, 'cdigit', p.cdigit, 'sourcetype', p.sourcetype, 'output', 'struct');
            im_photon = spcSMCphotonframe(smcstruct, 'combined', true, 'crop', smccrop);
    end
    
end

%% Make image
if dophotons
    % Reference
    if size(im_photon,3) > 1
        if p.usemean
            im = mean(im_photon,3);
        else
            im = median(double(im_photon),3);
        end
    else
        % im_photon is a sum
        im = double(im_photon) / n;
    end
    
    % Save a unfiltered copy
    if p.saveraw
        % Write
        writetiff(im, sprintf('%sraw.tif', fullfile(spcpaths.fp_out,spcpaths.cp_toseg(1:end-4))), 'double');
    end
    
    im(isnan(im)) = 0;
    
    % 2D Median filter
    if ~isempty(p.medfilt2size)
        % Apply median filter to reference
        im = medfilt2(im, p.medfilt2size, 'symmetric');
    end
    
    % Local normalize
    if p.uselocalnorm
        im_prime = single(im)-single(imgaussfilt(double(im), p.hp_norm_sigmas(1)));
        im = im_prime ./ (imgaussfilt(im_prime.^2, p.hp_norm_sigmas(2)) .^ (1/2));
        im(isnan(im)) = 0;
    end
    
    % Show
    if p.showresult
        figure;
        imshow(im, []);
    end
end 

%% Save
% Apply shifts to photon data
if dophotons
    % Write
    writetiff(im, fullfile(spcpaths.fp_out,spcpaths.cp_toseg), 'double');
end

end
