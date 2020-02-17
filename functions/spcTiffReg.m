function spcTiffReg(mouse, date, run, varargin)
% spcTiffReg register tiff of photon images, and apply the shifts to tm

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% File variables
addOptional(p, 'autogentiff', true); % If the original tiffs were not made, make them
addOptional(p, 'force', false); % Force overwrite or not

% Registration variables
addOptional(p, 'binxy', 1); % Binning

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

% Redo tm or not
if exist(fullfile(spcpaths.fp_out,spcpaths.regtif_tm), 'file')
    % File already exists
    if p.force
        % Force redo
        dotm = true;
    elseif input('Tm reg tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
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

% Make sure that the base photon tiff file is still around
if dophotons
    if ~exist(fullfile(spcpaths.fp_out,spcpaths.tif_photons), 'file')
        if p.autogentiff
            redotiffphotons = true;
        elseif input('No photon tiff file, make it? (1 = yes, 0 = no): ') == 1
            redotiffphotons = true;
        else
            redotiffphotons = false;
        end
    else
        redotiffphotons = false;
    end
else
    redotiffphotons = false;
end

% Make sure that the base tm tiff file is still around
if dotm
    if ~exist(fullfile(spcpaths.fp_out,spcpaths.tif_tm), 'file')
        if p.autogentiff
            redotifftm = true;
        elseif input('No tm tiff file, make it? (1 = yes, 0 = no): ') == 1
            redotifftm = true;
        else
            redotifftm = false;
        end
    else
        redotifftm = false;
    end

else
    redotifftm = false;
end

% Remake the base photonsa nd tm tiff files if needed
if redotiffphotons || redotifftm
    spcTiff(mouse, date, run, 'server', p.server, 'user', p.user, 'slice',...
        p.slice, 'cdigit', p.cdigit, 'photons', redotiffphotons, 'tm', redotifftm);
end

%% Register
if dophotons || dotm
    % Read
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_photons));
    if dotm
        im_tm = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_tm));
    end
    
    % Binning
    if p.binxy > 1
        im_photon_bin = binxy(im_photon, p.binxy);
    else
        im_photon_bin = im_photon;
    end
    
    % Median
    im2reg = median(im_photon_bin,3);
    
    % Get a focus area for calculating shifts
    imshow(im2reg,[]);
    h = imrect;
    regfocus = wait(h);
    regfocus = round(regfocus);
    close(gcf)
    
    % Initialize
    im_photon2reg = zeros(regfocus(4), regfocus(3), size(im_photon_bin,3));
    
    % Get coordinates
    regfocus2 = [regfocus(2), regfocus(2) + regfocus(4) - 1, regfocus(1), regfocus(1) + regfocus(3) - 1];
    
    % Crop reference
    im2reg = im2reg(regfocus2(1) : regfocus2(2), regfocus2(3) : regfocus2(4));
    
    % Normalize reference
    f_prime = im2reg - imgaussfilt(im2reg, 8);
    im2reg = f_prime ./ (imgaussfilt(f_prime.^2, 30).^(1/2));
    
    % Loop through and normalize/crop each frame
    for i = 1 : size(im_photon_bin, 3)
        frame = im_photon_bin(regfocus2(1) : regfocus2(2), regfocus2(3) : regfocus2(4), i);
        f_prime = frame - imgaussfilt(frame, 8);
        g_prime = f_prime ./ (imgaussfilt(f_prime.^2, 30).^(1/2));
        
        im_photon2reg(:,:,i) = g_prime;
    end
    
    % Get shifts
    [xy_shifts,~]=stackRegisterMA_RR(im_photon2reg,im2reg);
    
    % Modify the shifts according to binning (correlation values will be
    % messed up.
    if p.binxy > 1
        xy_shifts = xy_shifts * p.binxy;
    end
    
    % Apply shifts to photon data
    if dophotons
        % Apply shifts
        [~,im_photon]=stackRegisterMA_RR(im_photon, [], [], xy_shifts);
        
        % Write
        writetiff(im_photon, fullfile(spcpaths.fp_out,spcpaths.regtif_photons), 'double');
    
    end
    
    % Apply shifts to tm data
    if dotm
        % Apply shifts
        [~,im_tm]=stackRegisterMA_RR(im_tm, [], [], xy_shifts);
        
        % Write
        writetiff(im_tm, fullfile(spcpaths.fp_out,spcpaths.regtif_tm), 'double');
    end
end
end