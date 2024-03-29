function spcTiffDemonsReg(mouse, date, run, varargin)
% spcTiffDemonsReg uses demonsreg on the spc photon data and apply the 
% shifts to the tm data 

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

% Registration variables
addOptional(p, 'binxy', 1); % Binning
addOptional(p, 'iterations', 1); % Repeated registrations

% Edge
addOptional(p, 'edges', [0 0 0 0]); % Use edges to the processing (needed for bidirectional scanning).

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

% Redo tm or not
if exist(fullfile(spcpaths.fp_out,spcpaths.demregtif_tm), 'file')
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

% Redo tm or not
if exist(fullfile(spcpaths.fp_out,spcpaths.demregtif_iem), 'file')
    % File already exists
    if p.force
        % Force redo
        doiem = true;
    elseif input('IEM reg tiff file already exists, redo? (1 = yes, 0 = no): ') == 1
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

%% Load
if dophotons || dotm || doiem
    % Switch types
    switch p.sourcetype
        case 'raw'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.tif_photons);
            fp_tm = fullfile(spcpaths.fp_out, spcpaths.tif_tm);
            fp_iem = fullfile(spcpaths.fp_out, spcpaths.tif_iem);
        case 'registered'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.regtif_photons);
            fp_tm = fullfile(spcpaths.fp_out, spcpaths.regtif_tm);
            fp_iem = fullfile(spcpaths.fp_out, spcpaths.regtif_iem);
        case 'warped'
            fp_photon = fullfile(spcpaths.fp_out, spcpaths.warptif_photons);
            fp_tm = fullfile(spcpaths.fp_out, spcpaths.warptif_tm);
            fp_iem = fullfile(spcpaths.fp_out, spcpaths.warptif_iem);
    end
    
    % Read
    if exist(fp_photon, 'file')
        im_photon = readtiff(fp_photon);
    else
        disp('No photon tiff file. Check source type.');
        return;
    end
    if dotm
        if exist(fp_tm, 'file')
            im_tm = readtiff(fp_tm);
        else
            disp('No Tm input tiff file. Skip.');
            dotm = false;
        end
    end
    if doiem
        if exist(fp_iem, 'file')
            im_iem = readtiff(fp_iem);
        else
            disp('No IEM input tiff file. Skip.');
            dotm = false;
        end
    end
end

%% Iterations
for iter = 1 : p.iterations
    fprintf('Iterative registration: %i/%i\n', iter, p.iterations);
    
    %% Reference
    if dophotons || dotm
        % Reference
%         ref = median(im_photon,3);
        ref = mean(im_photon,3);
        
        % Apply edge
        if any(p.edges ~= 0)
            % Reference
            ref = ref(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :);

            if mod(size(ref,2), p.binxy) ~= 0 || mod(size(ref,1), p.binxy) ~= 0
                fprintf('Wrong edge values.\n')
                fprintf('Residual in the first two values: %i.\n', mod(size(ref,2),p.binxy));
                fprintf('Residual in the last two values: %i.\n', mod(size(ref,1), p.binxy));
                return;
            end

            % Photon data
            im_photon2 = im_photon(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :);
        else
            % Copy
            im_photon2 = im_photon;
        end
    end   

    %% Initialize
    if dophotons || dotm
        fprintf('Initializing arrays...')
        tic;
        
        if iter == 1 % Initialize first iteration only
            % Get the meshgrid
            [xx, yy] = meshgrid(1 : size(im_photon, 2), 1 : size(im_photon, 1));

            % Initialize registered data (using the grid file xx to get size)
            photon_reg = zeros(size(xx));
            photon_reg = repmat(photon_reg, [1, 1, size(im_photon,3)]);

            if dotm
                tm_reg = photon_reg;
            end
            if doiem
                iem_reg = photon_reg;
            end
        end
        
        % Apply median filter to reference
        ref = medfilt2(ref, p.medfilt2size, 'symmetric');

        % Bin ref, if necessary
        if p.binxy > 1
            ref = binxy(ref, p.binxy); 
        end

        % Highpass and normalize reference
        n = p.hp_norm_sigmas(1);
        m = p.hp_norm_sigmas(2);
        if p.uselocalnorm
            ref(isnan(ref)) = 0;
            ref_prime = single(ref)-single(imgaussfilt(double(ref),n));
            ref = ref_prime ./ (imgaussfilt(ref_prime.^2,m) .^ (1/2));
            ref(isnan(ref)) = 0;
        end
        figure;
        imshow(ref,[]);

        % Bin images
        if p.binxy > 1
            im_photon2 = binxy(im_photon2, p.binxy);
        end

        % Apply median filter and spatial normalization to images
        for  i = 1:size(im_photon2,3)
            % Prepare the data by median filtering and local normalizing
            if ~isempty(p.medfilt2size)
                im_photon2(:,:,i) = medfilt2(im_photon2(:,:,i), p.medfilt2size, 'symmetric');
            end
            if p.uselocalnorm
                f = im_photon2(:,:,i);
                f(isnan(f)) = 0;
                f_prime = f - imgaussfilt(single(f),n);
                g_prime = f_prime ./ (imgaussfilt(f_prime.^2,m).^(1/2));

                g_prime(isnan(g_prime)) = 0;

                im_photon2(:,:,i)=g_prime;
            end
        end

        % Initialize D_combined
        D_combined = zeros(size(im_photon2,1), size(im_photon2,2), size(im_photon2,3) * 2);
        t = toc;
        fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    end

    %% Register
    if dophotons || dotm
        fprintf('Registration...')
        tic;
        % Get the D tensor
        for  i = 1:size(im_photon2,3)
            % Demons reg
            [D,~] = imregdemons(im_photon2(:,:,i), ref, [32 16 8 4],...
                    'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitbar',false);

            D_combined(:, :, i*2-1 : i*2) = D; 
        end

        % resize back if necessary
        if p.binxy > 1
            D_combined = imresize(p.binxy * D_combined, p.binxy); % resize to bring back to full size of movie
        end

        % Re-embed D-combined into the full-size version if using edges
        if any(p.edges ~= 0)
            D_combined_full = zeros(size(im_photon,1), size(im_photon,2), size(im_photon,3) * 2);
            D_combined_full(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :) = D_combined;

            % Use the new D_combined for subsequent processing
            D_combined = D_combined_full;
            clear D_combined_full
        end

        % Data2 is a large variable
        clear im_photon2;
        t = toc;
        fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    end

    %% Apply the warp to original movie
    % Apply
    if dophotons || dotm
        % Reconstruct image stack
        fprintf('Reconstruct image stack...')
        for i = 1 : size(im_photon, 3)
            photon_reg(:,:,i) = ...
                interp2(xx, yy, im_photon(:,:,i),...
                xx + D_combined(:,:,2*i-1), yy + D_combined(:,:,2*i));

            if dotm
                tm_reg(:,:,i) = ...
                    interp2(xx, yy, im_tm(:,:,i),...
                    xx + D_combined(:,:,2*i-1), yy + D_combined(:,:,2*i));
            end
            
            if doiem
                iem_reg(:,:,i) = ...
                    interp2(xx, yy, im_iem(:,:,i),...
                    xx + D_combined(:,:,2*i-1), yy + D_combined(:,:,2*i));
            end
        end
        t = toc;
        fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    end
    
    %% Prep for next iteration
    if iter < p.iterations
        im_photon = photon_reg;
        
        if dotm
            im_tm = tm_reg;
        end
        if doiem
            im_iem = iem_reg;
        end
    end
end

%% Save
fprintf('Saving final stack...')
% Apply shifts to photon data
if dophotons
    % Write
    writetiff(photon_reg, fullfile(spcpaths.fp_out,spcpaths.demregtif_photons), 'double');
end

% Apply shifts to tm data
if dotm
    % Write
    writetiff(tm_reg, fullfile(spcpaths.fp_out,spcpaths.demregtif_tm), 'double');
end

% Apply shifts to iem data
if doiem
    % Write
    writetiff(iem_reg, fullfile(spcpaths.fp_out,spcpaths.demregtif_iem), 'double');
end
fprintf(' Done. \n');
   
end
