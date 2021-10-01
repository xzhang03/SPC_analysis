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
addOptional(p, 'uselocalnorm', true);
addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric); % Sigma for gaussian fit
addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter

% Show
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
        case 'registered'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.regtif_photons);
        case 'warped'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.warptif_photons);
        case 'demonsreg'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.demregtif_photons);
    end
    
    % Read
    im_photon = readtiff(fp_photon);
end

%% Make image
if dophotons
    % Reference
    im = median(im_photon,3);
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