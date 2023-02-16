function [tm, trace, photontrace] = spcGrossTm(mouse, date, run, varargin)
% spcGrossTm estimates gross Tm of a file
% [tm, trace, photontrace] = spcGrossTm(mouse, date, run, varargin)

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice

% File variables
addOptional(p, 'sourcetype', 'registered'); % Can be 'raw', 'reigstered', 'warped'

% Registration variables
addOptional(p, 'binxy', 1); % Binning

% Edge
addOptional(p, 'edges', [0 0 0 0]); % Use edges to the processing (needed for bidirectional scanning).
addOptional(p, 'ROI',false);

% Make plot
addOptional(p, 'makeplot', false); % Make histogram

% Drop out pixels
addOptional(p, 'dropoutthrehsold', 500); % Values below which are dropped out

% Smoothing for trace
addOptional(p, 'smoothmethod', @movmedian);
addOptional(p, 'smoothwin', 30);

% Photon trace
addOptional(p, 'photontrace', true);

% Figure position
addOptional(p, 'pos', [600 500 1000 420]);


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
    'slice', p.slice, 'cdigit', 1);

%% Load
% Switch types
switch p.sourcetype
    case 'raw'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.tif_photons); % in case use it later
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.tif_tm);
    case 'registered'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.regtif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.regtif_tm);
    case 'warped'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.warptif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.warptif_tm);
    case 'demonsreg'
        fp_photon = fullfile(spcpaths.fp_out, spcpaths.demregtif_photons);
        fp_tm = fullfile(spcpaths.fp_out, spcpaths.demregtif_tm);
end

% Read
im_tm = readtiff(fp_tm);
if p.photontrace || p.ROI
    im_photon = readtiff(fp_photon);
end

%% Preprocessing
% Edge
if any(p.edges ~= 0)
    im_tm = im_tm(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :);
end

% Bin
im_tm = binxy(im_tm, p.binxy);
if p.photontrace || p.ROI
    im_photon = binxy(im_photon, p.binxy);
end
%% ROI
if p.ROI
    % Read
    ROI = getpoly(mean(im_photon,3));
    for i = 1 : size(im_tm,3)
        f = im_tm(:,:,i);
        f(~ROI) = nan;
        im_tm(:,:,i) = f;
        
        if p.photontrace
            f = im_photon(:,:,i);
            f(~ROI) = nan;
            im_photon(:,:,i) = f;
        end
    end
end

%% Reshape
% Reshape
[x, y, z] = size(im_tm);
im_tm2 = reshape(im_tm, [x*y*z 1]);
im_tm2 = im_tm2(im_tm2 >= p.dropoutthrehsold);

tm = nanmean(im_tm2);

%% Get timecourse
% Time course
trace = zeros(z,1);
for i = 1 : z
    f = im_tm(:,:,i);
    trace(i) = nanmean(f(f >= p.dropoutthrehsold));
end
if p.smoothwin > 1
    trace = p.smoothmethod(trace, p.smoothwin);
end

%% Get photon timecourse
if p.photontrace
    % Time course
    photontrace = zeros(z,1);
    for i = 1 : z
        photontrace(i) = nanmean(nanmean(im_photon(:,:,i)));
    end
    if p.smoothwin > 1
        photontrace = p.smoothmethod(photontrace, p.smoothwin);
    end
else
    photontrace = [];
end

%% Plot
if p.makeplot
    % Panels
    npanels = 2+p.photontrace;
    
    % Histogram
    figure('Position', p.pos)
    subplot(1,npanels,1)
    hist(im_tm2, 1000);
    title(sprintf('%s %s run%i Tm histogram', mouse, date, run));
    
    % Show mean
    hold on
    ylims = get(gca, 'ylim');
    plot([tm tm], ylims, 'r-', 'LineWidth', 3)
    hold off
    
    text(double(tm)*1.02, ylims(2)*0.97, sprintf('Tm = %i', round(tm)), ...
        'FontSize', 15, 'Color', 'r');
    
    % Subplot 2
    subplot(1,npanels,2);
    plot(trace)
    
    if p.photontrace
        subplot(1,npanels,3);
        plot(photontrace)
    end
    
end

end
