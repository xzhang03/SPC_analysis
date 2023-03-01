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
addOptional(p, 'autoc', false); % Automate cdigit - experimental
addOptional(p, 'frommat', true); % Loda from mat if exists

% Imaage unscrambling
addOptional(p, 'rowinc', 2044.72); % This is the key variable to descramble images. Fine adjust for shearing

% Time domain: tm and iem
addOptional(p, 'tbins', 256); % Time bins
addOptional(p, 'tcycle', 12500); % in ps
addOptional(p, 'T1', 21); % data before this is not considered for tm and iem
addOptional(p, 'T2', 240); % data after this is not considered for tm and iem

% Resize (non-linear, correct for warping)
addOptional(p, 'resize', true);
addOptional(p, 'resizedim', [511 1250]);

% IRF
addOptional(p, 'deconvforiem', false); % If set true it takes 10x as long
addOptional(p, 'irf', [882; 6176; 15000; 15000; 3529]);

% Output
addOptional(p, 'save', true);
addOptional(p, 'makeplot', true);
addOptional(p, 'pos', [600 500 1000 420]);
addOptional(p, 'smoothwin', 3); % Change to 0 for no smoothing

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

%% Initialize
if p.frommat
    if exist(fullfile(spcpaths.fp_out, spcpaths.peek), 'file')
        savestruct = load(fullfile(spcpaths.fp_out, spcpaths.peek));
        fprintf('Loaded from previous mat.\n');
        donew = false;
        p.save = false;
    else
        donew = true;
    end
else
    donew = true;
end

if donew
    % Traces
    photontrace = zeros(spcpaths.n, 1);
    tmtrace = zeros(spcpaths.n, 1);
    iemtrace = zeros(spcpaths.n, 1);

    % Decays
    tracemat = zeros(p.tbins, spcpaths.n);
    tracedcmat = zeros(p.tbins, spcpaths.n);
else
    % Traces
    photontrace = savestruct.photontrace;
    tmtrace = savestruct.tmtrace;
    iemtrace = savestruct.iemtrace;

    % Decays
    tracemat = savestruct.tracemat;
    tracedcmat = savestruct.tracedcmat;
end

%% Load
if donew
    hwait = waitbar(0);
    for ind = spcpaths.cinds
        waitbar(ind/spcpaths.n, hwait, sprintf('Processing %s %s run%i: %i/%i', mouse, date, run, ind, spcpaths.n));

        % Load movie
        tic
        if p.resize
            mov = spcLoadsdt(fullfile(spcpaths.fp, sprintf(spcpaths.sdt_in,ind)), p.rowinc, p.resizedim);
        else
            mov = spcLoadsdt(fullfile(spcpaths.fp, sprintf(spcpaths.sdt_in,ind)), p.rowinc, []);
        end
        t = toc;
        fprintf('Loaded Frame %i in %0.2f s.\n', ind, t);

        % Trace
        trace = squeeze(sum(sum(mov,1),2));
        trace = trace(ivec);

        % Photon count
        photontrace(ind) = sum(trace);

        % Tm
        tmtrace(ind) = sum(trace .* tvec) / sum(trace);

        % IEM
        tracedc = deconvlucy(trace,p.irf);
        iemtrace(ind) = sum(tracedc) / max(tracedc) * tres;

        % Save traces
        tracemat(:,ind) = cat(1, zeros(p.T1-1, 1), trace, zeros((p.tbins-p.T2), 1));
        tracedcmat(:,ind) = cat(1, zeros(p.T1-1, 1), tracedc, zeros((p.tbins-p.T2), 1));
    end
    close(hwait)
end

%% Plot
if p.makeplot
    % Panels
    npanels = 3;
    
    % Photon
    figure('Position', p.pos)
    subplot(1,npanels,1)
    if p.smoothwin > 0
        plot(smooth(photontrace, p.smoothwin));
    else
        plot(photontrace);
    end
    title(sprintf('%s %s run%i Photons', mouse, date, run));
    
    % Tm
    subplot(1,npanels,2)
    if p.smoothwin > 0
        plot(smooth(tmtrace, p.smoothwin));
    else
        plot(tmtrace);
    end
    title(sprintf('%s %s run%i Tm', mouse, date, run));
    
    % IEM
    subplot(1,npanels,3)
    if p.smoothwin > 0
        plot(smooth(iemtrace, p.smoothwin));
    else
       plot(iemtrace);
    end
    title(sprintf('%s %s run%i IEM', mouse, date, run));
end

%% Save
if p.save
    % Save struct
    savestruct = struct('photontrace', photontrace, 'tmtrace', tmtrace, 'iemtrace', iemtrace, 'tracemat',...
        tracemat, 'tracedcmat', tracedcmat, 'p', p);
    
    save(fullfile(spcpaths.fp_out, spcpaths.peek), '-struct', 'savestruct', '-v7.3');
    disp('Saved');
end

end