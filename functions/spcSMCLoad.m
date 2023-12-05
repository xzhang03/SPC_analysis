function [output, n, smccrop] = spcSMCLoad(mouse, date, run, varargin)
% Load smc structures

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)
addOptional(p, 'autoc', true); % Auto c digit
addOptional(p, 'spcpaths', {}); % You can just pass the spcpaths struct here

% File variables
addOptional(p, 'sourcetype', 'smc'); % Can be 'smc', 'smcreg', 'smcdereg'
addOptional(p, 'output', 'struct'); % Can be 'cell' or 'struct'
addOptional(p, 'usewaitbar', true);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
if isempty(p.spcpaths)
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

    % IO
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
else
    spcpaths = p.spcpaths;
end
n = spcpaths.n;

%% Reading
% Initialize cell
tic;
fprintf('Loading %i %s... ', n, p.sourcetype);
smccell = cell(n,1);
if p.usewaitbar
    hwaitsmc = waitbar(0, sprintf('Loading %s %i/%i', p.sourcetype, 1, n));
end
for i = 1 : n
    if mod(i,5) == 0 && p.usewaitbar
        waitbar(i/n, hwaitsmc, sprintf('Loading %s %i/%i', p.sourcetype, i, n));
    end
    loaded = load(fullfile(spcpaths.fp, sprintf(spcpaths.(p.sourcetype), i)), '-mat');
    smccell{i} = loaded.smc;
    if i == 1
        if isfield(loaded, 'crop')
            smccrop = loaded.crop;
        else
            smccrop = [];
        end
    end
end
if p.usewaitbar
    close(hwaitsmc);
end
t = toc;
fprintf('Done. %0.1f s.\n', t);

if strcmpi(p.output, 'cell')
    output = smccell;
end

%% Convert to struct if necessary
% Convert to struct
if strcmpi(p.output, 'struct')
    l = length(smccell{1});
    smstruct = cat(1, smccell{:});
        
    % Counter
    counter = 0;
    ind2 = 1;
    for i = 1 : length(smstruct)
        counter = counter + 1;
        if counter > l
            counter = 1;
            ind2 = ind2 + 1;
        end
        smstruct(i).ind2 = ind2;
    end
    output = smstruct;
end
end