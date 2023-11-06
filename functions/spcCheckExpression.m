function spcCheckExpression(mouse, date, run, varargin)
%spcCheckBinsmc Checks xy bins for lifetime calculations

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Parameters
addOptional(p, 'compress', true); % First time it will make compression file, later reading it
addOptional(p, 'frames', 1:5);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;


%% IO
% Get paths
spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);

%% Load
hwait = waitbar(0);
for ind = p.frames
    waitbar(ind/p.frames, hwait, sprintf('Processing %s %s run%i: %i/%i', mouse, date, run, ind, p.frames));

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

    % Initialize
    if ind == 1
        % Get size
        sz = size(mov);
    
        % Initialize data
        fprintf('Intializing... ');
        tic;
        mov3d = uint16(zeros(sz(1), sz(2), length(p.frames)));
        d = whos('mov3d');
        t = toc;
        fprintf('%0.2f GB in RAM in %0.2f s.\n', d.bytes / 1e9, t);
    end

    % Trim and save
    mov3d(:,:,ind) = uint16(sum(mov,3));

    t = toc;
    fprintf(' Processing %0.2f s.\n', t);

end
close(hwait)

%% Plot
figure();
imshow(mean(mov3d,3), []);

end

