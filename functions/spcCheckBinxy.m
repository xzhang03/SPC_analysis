function spcCheckBinxy(mouse, date, run, varargin)
%spcCheckBinsmc Checks xy bins for lifetime calculations

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Parameters
addOptional(p, 'bins2check', 1:2:11);
addOptional(p, 'frame2check', 2);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

% N Bins to check
nchecks = length(p.bins2check);

%% IO
% Get paths
spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit);
if exist(fullfile(spcpaths.fp, sprintf(spcpaths.smc, p.frame2check)), 'file')
    smc = load(fullfile(spcpaths.fp, sprintf(spcpaths.smc, p.frame2check)), '-mat');
    im = spcDecompress(smc.smc);
else
    disp('SMC not found.')
    return;
end

%% Get input
% Show and get input
im2show = sum(im,3);
figure;
imshow(im2show, []);
[x, y] = ginput(1);
close(gcf);

% Round
x = round(x);
y = round(y);

%% Plot
figure('Position', [50,50,1600,500]);
for i = 1 : nchecks
    subplot(1, nchecks, i);
    t = spcMovtrace(im, y, x, p.bins2check(i));
    plot(t);
    title(sprintf('XYBIN = %i', p.bins2check(i)));
end

end

