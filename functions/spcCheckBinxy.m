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
addOptional(p, 'binxy2check', 1:2); % Frame xy bin
addOptional(p, 'bins2check', 1:2:11); % Lifetime xy bin
addOptional(p, 'frame2check', 2);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

% N Bins to check
nframebins = length(p.binxy2check);
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
for ii = 1 : nframebins
    if p.binxy2check(ii) == 1
        im2 = im;
        x2 = x;
        y2 = y;
    else
        im2 = binxy(im, p.binxy2check(ii)) * (p.binxy2check(ii))^2;
        x2 = round(x / p.binxy2check(ii));
        y2 = round(y / p.binxy2check(ii));
    end
    for i = 1 : nchecks
        ind = (p.binxy2check(ii)-1) * nchecks + i;
        subplot(nframebins, nchecks, ind);
        t = spcMovtrace(im2, y2, x2, p.bins2check(i));
        plot(t);
        title(sprintf('XY Bin = %i | Local Bin = %i', p.binxy2check(ii), p.bins2check(i)));
    end
end

end

