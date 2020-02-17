function spcTiff(mouse, date, run, varargin)
% spcTiff converts spc outputs to tiff files

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% File variables
addOptional(p, 'photons', true); % Process photons file or not
addOptional(p, 'tm', true); % Process tm file or not
addOptional(p, 'force', false); % Force overwrite or not

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

% Do photon data?
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

%% Photons
if dophotons
    % Read one frame and initialize
    ind = spcpaths.cinds(1);
    im_photon = load(fullfile(spcpaths.fp, sprintf(spcpaths.photons_in,ind)));
    im_photon = repmat(im_photon, [1 1 length(spcpaths.cinds)]);
    
    for i = 2 : length(spcpaths.cinds)
        % Filename
        ind = spcpaths.cinds(i);
        fn_photon = fullfile(spcpaths.fp, sprintf(spcpaths.photons_in,ind));
        
        % Load
        im_photon(:,:,ind) = load(fn_photon);
    end
    
    % Write
    writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.tif_photons), 'double');
end

%% Tm
if dotm
    % Read one frame and initialize
    ind = spcpaths.cinds(1);
    im_tm = load(fullfile(spcpaths.fp, sprintf(spcpaths.tm_in,ind)));
    im_tm = repmat(im_tm, [1 1 length(spcpaths.cinds)]);
    
    for i = 2 : length(spcpaths.cinds)
        % Filename
        ind = spcpaths.cinds(i);
        fn_tm = fullfile(spcpaths.fp, sprintf(spcpaths.tm_in,ind));
        
        % Load
        im_tm(:,:,ind) = load(fn_tm);
    end
    
    % Write
    writetiff(im_tm, fullfile(spcpaths.fp_out,spcpaths.tif_tm), 'double');
end
end