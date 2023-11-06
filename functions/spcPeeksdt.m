function savestruct = spcPeeksdt(mouse, date, run, varargin)
% spcPeeksdt superficially peeks at raw sdt data
% savestruct = spcPeeksdt(mouse, date, run, varargin)

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)
addOptional(p, 'autoc', true); % Automate cdigit - experimental
addOptional(p, 'frommat', true); % Loda from mat if exists

% Compress
addOptional(p, 'autosave', false);
addOptional(p, 'compress', true); % First time it will make compression file, later reading it

% Time domain: tm and iem
addOptional(p, 'tbins', 256); % Time bins
addOptional(p, 'tcycle', 12500); % in ps
addOptional(p, 'T1', 21); % data before this is not considered for tm and iem
addOptional(p, 'T2', 240); % data after this is not considered for tm and iem

% Cropping
addOptional(p, 'usecrop', true);
addOptional(p, 'autocrop', false);
addOptional(p, 'autocrop_size', 3);
addOptional(p, 'autocrop_threshold', 1);
addOptional(p, 'showautocrop', true); % Show autocrop if using it
addOptional(p, 'crop', []);

% Binning (sacrificing spatial resolution and noise separability)
addOptional(p, 'binxy', 1);
addOptional(p, 'bint', 1);

% Default thresholding
addOptional(p, 'thresh', [0 9000]); % Default tm thresh
addOptional(p, 'phthresh', [0 1000]);  % Default photonthresh

% IRF
addOptional(p, 'deconvforiem', false); % If set true it takes 10x as long
addOptional(p, 'irf', [882; 6176; 15000; 15000; 3529]);

% Output
addOptional(p, 'makeplot', true);
addOptional(p, 'pos', [50 500 1800 420]);
addOptional(p, 'smoothwin', 3); % Change to 0 for no smoothing
addOptional(p, 'tmLUT', [0 3500]);

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

fovtrace = [];
ROItrace = [];

%% Time vectors
% Time resolution
tres = p.tcycle / p.tbins * p.bint;

% Time vector
ivec = round(p.T1/p.bint : p.T2/p.bint)';
tvec = (ivec - ivec(1)) * tres;

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

%% Initialize
if p.frommat
    if exist(fullfile(spcpaths.fp_out, spcpaths.peek), 'file')
        % If not making plots stop here
        if ~p.makeplot
            fprintf('Not making new structs or new plots. Stopping here.\n');
            return;
        end
        
        savestruct = load(fullfile(spcpaths.fp_out, spcpaths.peek));
        fprintf('Loaded from previous mat.\n');
        donew = false;
        
    else
        donew = true;
    end
else
    donew = true;
end

if ~donew
    % Save changeable variables
    makeplot = p.makeplot;
    fprintf('Decompressing... ');
    tic;
    mov4d = spcDecompress(savestruct.mov_compressed);
    t = toc;
    fprintf('Done. %0.2f s\n', t);
    tm3d = savestruct.tm3d;
    tres = savestruct.tres;
    ivec = savestruct.ivec;
    tvec = savestruct.tvec;
    mouse = savestruct.mouse;
    date = savestruct.date;
    run = savestruct.run;
    p = savestruct.p;
    
    % Load changeable variables
    p.makeplot = makeplot;
    
    % Not autosave if loaded
    p.autosave = false;
end

%% Catch errors
if isempty(spcpaths.cinds)
    fprintf('Found no sdt files in %s\n', spcpaths.fp);
    if exist(spcpaths.fp, 'dir')
        fprintf('Path exists.\n');
    else
        fprintf('Path does not exists.\n');
    end
    return;
end

%% Load
if donew
    hwait = waitbar(0);
    for ind = spcpaths.cinds
        waitbar(ind/spcpaths.n, hwait, sprintf('Processing %s %s run%i: %i/%i', mouse, date, run, ind, spcpaths.n));
        
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
            mov = spcDecompress(smc.smc);
        else
            mov = spcLoadsdt(fullfile(spcpaths.fp, sprintf(spcpaths.sdt_in,ind)));
            
            if p.compress
                % Save smc for future
                smc = spcCompress(mov);
                save(fullfile(spcpaths.fp, sprintf(spcpaths.smc,ind)), 'smc', '-v7.3');
            end
        end
        t = toc;
        fprintf('Loaded Frame %i in %0.2f s.', ind, t);
        
        % Cropping
        if p.usecrop && ind == 1
            if isempty(p.crop)
                if p.autocrop
                    p.crop = spcAutocrop(sum(mov,3), 'size', p.autocrop_size, 'threshold', p.autocrop_threshold);
                    if p.showautocrop
                        figure
                        imshow(sum(mov,3),[]);
                        title(sprintf('%s %s run%i Autocrop', mouse, date, run))
                        rectangle('Position', [p.crop(1) p.crop(2) p.crop(3)-p.crop(1) p.crop(4)-p.crop(2)],...
                            'EdgeColor','g', 'LineWidth',2);
                    end
                else
                    figure
                    imshow(sum(mov,3),[]);
                    p.crop = round(wait(imrect()));
                    p.crop(3) = min(p.crop(1) + p.crop(3), size(mov,2));
                    p.crop(1) = max(p.crop(1),1);
                    p.crop(4) = min(p.crop(2) + p.crop(4), size(mov,1));
                    p.crop(2) = max(p.crop(2),1);
                    close(gcf)

                    disp(p.crop);
                end
            end
        end
        mov = mov(p.crop(2):p.crop(4), p.crop(1):p.crop(3), :);
        
        tic
        % Binning
        if p.binxy > 1
            mov = binxy(mov, p.binxy);
            mov = mov * p.binxy * p.binxy;
        end
        if p.bint > 1
            mov = bint(mov, p.bint);
            mov = mov * p.bint;
        end
        
        % Get size
        sz = size(mov);
        
        % Initialize
        if ind == 1
            % Initialize data
            fprintf('Intializing... ');
            tic;
            mov4d = uint16(zeros(sz(1), sz(2), length(ivec), length(spcpaths.cinds)));
            d = whos('mov4d');
            t = toc;
            fprintf('%0.2f GB in RAM in %0.2f s.\n', d.bytes / 1e9, t);
            
            % Initialize tm
            ttensor = reshape(tvec, [1 1 length(tvec)]);
            ttensor = repmat(ttensor, [sz(1) sz(2)]);
            tm3d = zeros(sz(1), sz(2), length(spcpaths.cinds));
        end
        
        % Trim and save
        movtrim = mov(:, :, ivec);
        mov4d(:,:,:,ind) = uint16(movtrim);
        
        % Calculate field-wide tm
        movtrim = double(movtrim);
        tm3d(:,:,ind) = sum(movtrim .* ttensor, 3) ./ sum(movtrim, 3);
        
        t = toc;
        fprintf(' Processing %0.2f s.\n', t);
        
    end
    close(hwait)
end

%% Prepare for plotting
% Threshold
[mov4dthresh, tm3dthresh, failmask] = applytmthresh(mov4d, tm3d, p.thresh, p.phthresh);

% Calculations
[fov,tmfov, photontrace, tmtrace, iemtrace, decay] = datacal(mov4dthresh, tm3dthresh, p.irf, tvec, tres);

% RGB fov
rgbfov = makergbfov(fov, failmask);

%% Autosave (batching)
if p.autosave
    savedata();
end

%% Plot
if ~p.makeplot
    return;
end

% Panels
npanels = 8;
hfig = figure('Position', p.pos);

% Photon image
subplot(1,npanels,1:2)
hfov = imshow(rgbfov);
title(sprintf('%s %s run%i Photons', mouse, date, run));

% tm image
subplot(1,npanels,3:4)
htmfov = imshow(tmfov, []);
if isempty(p.tmLUT)
    p.tmLUT = round(htmfov.Parent.CLim);
end
title(sprintf('%s %s run%i Tm', mouse, date, run));

% Photon
subplot(1,npanels,5)
if p.smoothwin > 0
    hph = plot(smooth(photontrace, p.smoothwin));
else
    hph = plot(photontrace);
end
title(sprintf('Photons'));

% Tm
subplot(1,npanels,6)
if p.smoothwin > 0
    htm = plot(smooth(tmtrace, p.smoothwin));
else
    htm = plot(tmtrace);
end
title(sprintf('Tm: %0.0f', mean(tmtrace)));

% IEM
subplot(1,npanels,7)
if p.smoothwin > 0
   hiem = plot(smooth(iemtrace, p.smoothwin));
else
   hiem = plot(iemtrace);
end
title(sprintf('IEM: %0.0f', mean(iemtrace)));

% Decay
subplot(1,npanels,8)
hdecay = plot(tvec(1:40), decay(1:40));
title(sprintf('Decay'));

% Controls
% Panel
hpan = uipanel(hfig,'Position',[0.0100 0.1000 0.0900 0.8500]);

% Tm LUT
vp1 = [20 320 100 20];
uicontrol(hpan, 'Style', 'text', 'Position', vp1, 'String', 'Tm LUT:')
uicontrol(hpan, 'Style', 'text', 'Position', vp1 - [0 22 0 0], 'String', '-')
htmLUT1 = uicontrol(hpan,'Style','edit','Position', vp1 - [0 20 60 0], 'String', num2str(p.tmLUT(1)),...
    'callback', @tmLUTchange);
htmLUT2 = uicontrol(hpan,'Style','edit','Position', vp1 - [-60 20 60 0], 'String', num2str(p.tmLUT(2)),...
    'callback', @tmLUTchange);

% Photon thresh
vp2 = [20 260 100 20];
uicontrol(hpan, 'Style', 'text', 'Position', vp2, 'String', 'Photon thresh:')
uicontrol(hpan, 'Style', 'text', 'Position', vp2 - [0 22 0 0], 'String', '-')
hpthresh1 = uicontrol(hpan,'Style','edit','Position', vp2 - [0 20 60 0], 'String', num2str(p.phthresh(1)),...
    'callback', @threshchange);
hpthresh2 = uicontrol(hpan,'Style','edit','Position', vp2 - [-60 20 60 0], 'String', num2str(p.phthresh(2)),...
    'callback', @threshchange);


% Tm thresh
vp3 = [20 200 100 20];
uicontrol(hpan, 'Style', 'text', 'Position', vp3, 'String', 'Tm thresh:')
uicontrol(hpan, 'Style', 'text', 'Position', vp3 - [0 22 0 0], 'String', '-')
hthresh1 = uicontrol(hpan,'Style','edit','Position', vp3 - [0 20 60 0], 'String', num2str(p.thresh(1)),...
    'callback', @threshchange);
hthresh2 = uicontrol(hpan,'Style','edit','Position', vp3 - [-60 20 60 0], 'String', num2str(p.thresh(2)),...
    'callback', @threshchange);

% Status
vp4 = [20 140 100 20];
hwaitmsg = uicontrol(hpan, 'Style', 'text', 'Position', vp4, 'String', '', 'FontSize', 12);

% ROI
vp5 = [20 90 100 20];
uicontrol(hpan,'Style','pushbutton','Position', vp5, 'String',...
    'ROI', 'callback', @checkroi);

% Save screenshot
uicontrol(hpan,'Style','pushbutton','Position', vp5 - [0 40 0 0], 'String',...
    'Screenshot', 'callback', @savescreenshot);

% Save
uicontrol(hpan,'Style','pushbutton','Position', vp5 - [0 60 0 0], 'String',...
    'Save', 'callback', @savedata);

% Quit
uicontrol(hpan,'Style','pushbutton','Position', vp5 - [0 80 0 0], 'String',...
    'Quit', 'callback', @quit);

%% Anonymous functions
    % Anonymous functions
    function checkroi(~,~)
        figure(hfig);
        subplot(1,npanels,1:2)
        hrect = imrect();
        coords = round(wait(hrect));
        coords(3) = min(coords(1) + coords(3), size(fov,2));
        coords(1) = max(coords(1),1);
        coords(4) = min(coords(2) + coords(4), size(fov,1));
        coords(2) = max(coords(2),1);
        delete(hrect);
        
        % New figure
        figure('Position', p.pos - [0 p.pos(4) 0 0]);
        
        % Fov new
        subplot(1,npanels,1:2)
        imshow(rgbfov);
        rectangle('Position', [coords(1) coords(2) coords(3)-coords(1) coords(4)-coords(2)],...
            'EdgeColor','g', 'LineWidth',2)
        title(sprintf('%s %s run%i Photons', mouse, date, run));
        
        % tm image
        subplot(1,npanels,3:4)
        htmfov = imshow(tmfov, p.tmLUT);
        rectangle('Position', [coords(1) coords(2) coords(3)-coords(1) coords(4)-coords(2)],...
            'EdgeColor','g', 'LineWidth',2)
        title(sprintf('%s %s run%i Tm', mouse, date, run));
        
        % Cropping
        mov4d_crop = mov4d(coords(2):coords(4), coords(1):coords(3), :, :);
        mov4dthresh_crop = mov4dthresh(coords(2):coords(4), coords(1):coords(3), :, :);
        tm3d_crop = tm3d(coords(2):coords(4), coords(1):coords(3), :);
        tm3dthresh_crop = tm3dthresh(coords(2):coords(4), coords(1):coords(3), :);
        
        % Calculate cropped
        [~,~, photontrace_crop, tmtrace_crop, iemtrace_crop, decay_crop] = ...
            datacal(mov4d_crop, tm3d_crop, p.irf, tvec, tres);
        [~,~, photontrace_thresh_crop, tmtrace_thresh_crop, iemtrace_thresh_crop, decay_thresh_crop] = ...
            datacal(mov4dthresh_crop, tm3dthresh_crop, p.irf, tvec, tres);
        
        % Photon
        subplot(1,npanels,5)
        if p.smoothwin > 0
            plot(movmean([photontrace_thresh_crop, photontrace_crop], p.smoothwin, 1));
        else
            plot([photontrace_thresh_crop, photontrace_crop]);
        end
        title(sprintf('Photons'));

        % Tm
        subplot(1,npanels,6)
        if p.smoothwin > 0
            plot(movmean([tmtrace_thresh_crop', tmtrace_crop'], p.smoothwin, 1));
        else
            plot([tmtrace_thresh_crop', tmtrace_crop']);
        end
        title(sprintf('Tm'));

        % IEM
        subplot(1,npanels,7)
        if p.smoothwin > 0
           plot(movmean([iemtrace_thresh_crop, iemtrace_crop], p.smoothwin, 1));
        else
           plot([iemtrace_thresh_crop, iemtrace_crop]);
        end
        title(sprintf('IEM'));
        
        % ROI traces
        ROItrace = struct('coords', coords, 'photontrace_thresh_crop', photontrace_thresh_crop,...
            'tmtrace_thresh_crop', tmtrace_thresh_crop, 'iemtrace_thresh_crop', iemtrace_thresh_crop);

        % Decay
        subplot(1,npanels,8)
        plot(tvec(1:40), [decay_thresh_crop(1:40), decay_crop(1:40)]);
        title(sprintf('Decay'));
        legend({'Thresholded', 'Raw'}, 'Location', 'northeast');
    end
    
    % Apply threshold
    function [mov4d, tm3d, failmask] = applytmthresh(mov4d, tm3d, tmthresh, phthresh)
        fprintf('Thresholding... ');
        tic;
        
        % Get what passes
        tm3dfail1 = (tm3d < tmthresh(1)) | (tm3d > tmthresh(2));
        fovhere = sum(sum(mov4d,3), 4);
        tm3dfail2 = (fovhere < phthresh(1)) | (fovhere > phthresh(2));
        tm3dfail = tm3dfail1 | tm3dfail2;
        
        % Reshape passes
        mov4dfail = reshape(tm3dfail, [size(tm3dfail,1), size(tm3dfail,2), 1, size(tm3dfail,3)]);
        mov4dfail = repmat(mov4dfail, [1 1 size(mov4d, 3) 1]);

        % Apply to tm3d
        tm3d(tm3dfail) = nan;
        mov4d(mov4dfail) = 0;
        failmask = mean(tm3dfail, 3);
        
        t = toc;
        fprintf('Done. %0.2f s\n', t);
    end

    % Main calculation
    function [fov,tmfov, photontrace, tmtrace, iemtrace, decay] = datacal(mov4, tm3d, irf, tvec, tres)
        fprintf('Calculating... ');
        tic;
        
        % Get fov
        fov = sum(sum(mov4,3), 4);

        % Get tm fov
        tmfov = nanmean(tm3d, 3);

        % Mov 2d
        mov2d = squeeze(sum(sum(mov4, 1),2));
        
        % Decay
        decay = mean(mov2d, 2);
        
        % Tvec 2d
        tvec2d = tvec * ones(1, size(mov2d,2));

        % Photon trace
        photontrace = squeeze(sum(sum(sum(mov4,1),2),3));

        % Tm traces
        tmtrace = sum(mov2d .* tvec2d, 1) ./ sum(mov2d,1);

        % Mov 2d deconvolve
        iemtrace = zeros(size(mov2d,2),1);
        for i = 1 : size(mov2d,2)
            tracedc = deconvlucy(mov2d(:,i), irf);
            iemtrace(i) = sum(tracedc) / max(tracedc) * tres;
        end
        
        % trace struct
        fovtrace = struct('photon', photontrace, 'tm', tmtrace, 'iem', iemtrace);
        
        t = toc;
        fprintf('Done. %0.2f s\n', t);
    end

    % make RGB fov
    function rgbfov = makergbfov(fov, failmask)
        fov = mat2gray(fov);
        rgbfov = repmat(fov, [1 1 3]);
        rgbfov(:,:,1) = failmask;
    end

    % Change threshold
    function threshchange(~,~)
        hwaitmsg.String = 'Calculating...';
        drawnow();
        
        % Change threshold values
        thresh1 = str2double(hthresh1.String);
        thresh2 = str2double(hthresh2.String);
        p.thresh = [thresh1, thresh2];
        
        % Change photon threshold values
        phthresh1 = str2double(hpthresh1.String);
        phthresh2 = str2double(hpthresh2.String);
        p.phthresh = [phthresh1, phthresh2];
        
        % Threshold
        [mov4dthresh, tm3dthresh, failmask] = applytmthresh(mov4d, tm3d, p.thresh, p.phthresh);

        % Calculations
        [fov,tmfov, photontrace, tmtrace, iemtrace, decay] = datacal(mov4dthresh, tm3dthresh, p.irf, tvec, tres);
        
        % RGB fov
        rgbfov = makergbfov(fov, failmask);

        % Update
        hfov.CData = rgbfov;
        htmfov.CData = tmfov;
        
        % Update traces
        if p.smoothwin > 0
            photontrace = smooth(photontrace, p.smoothwin);
        end
        hph.YData = photontrace;
        hph.Parent.YLim = [min(photontrace) max(photontrace)];
        
        if p.smoothwin > 0
            tmtrace = smooth(tmtrace, p.smoothwin);
        end
        htm.YData = tmtrace;
        htm.Parent.YLim = [min(tmtrace) max(tmtrace)];
        
        if p.smoothwin > 0
            iemtrace = smooth(iemtrace, p.smoothwin);
        end
        hiem.YData = iemtrace;
        hiem.Parent.YLim = [min(iemtrace) max(iemtrace)];
        
        hdecay.YData = decay(1:40);
        hdecay.Parent.YLim = [min(decay(1:40)) max(decay(1:40))];
        
        hwaitmsg.String = '';
    end

    % Photon LUT
    function tmLUTchange(~, ~)
        LUTlower = str2double(htmLUT1.String);
        LUTupper = str2double(htmLUT2.String);
        p.tmLUT = [LUTlower LUTupper];
        htmfov.Parent.CLim = p.tmLUT;
    end
    
    function savescreenshot(~,~)
        [~, peekname, ~] = fileparts(spcpaths.peek);
        peekname_full = fullfile(spcpaths.fp_out, sprintf('%s.png', peekname));
        saveas(hfig, peekname_full, 'png');
    end

    function savedata(src, ~)
        hwaitmsg.String = 'Saving...';
        drawnow();

        % Save struct
        savestruct = struct('mov_compressed', spcCompress(mov4d), 'tm3d', tm3d, 'tres', tres, 'ivec', ivec,...
            'tvec', tvec, 'mouse', mouse, 'date', date, 'run', run, 'p', p, 'fovtrace', fovtrace, 'ROItrace', ROItrace);

        save(fullfile(spcpaths.fp_out, spcpaths.peek), '-struct', 'savestruct', '-v7.3');
        disp('Saved');

        hwaitmsg.String = '';
        src.String = 'Data saved';
    end

    function quit(~, ~)
        close(hfig);
    end
end