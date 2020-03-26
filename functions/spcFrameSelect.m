function spcFrameSelect(mouse,date,run, varargin)
% spcFrameSelect uses a UI to select which FLIM frames to include, so to
% remove the blury ones.

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Use registered?
addOptional(p, 'useregistered', true);

% Number of panels per row
addOptional(p, 'panelsperrow', 3);

% Binning
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

% Load
if p.useregistered
    % Load
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.regtif_photons));
    
    % Figure name
    fign = spcpaths.regtif_photons;
else
    % Load
    im_photon = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_photons));
    
    % Figure name
    fign = spcpaths.tif_photons;
end

% Number of frames
nframes = size(im_photon, 3);

% Binning if needed
if p.binxy > 1
    im_photon = binxy(im_photon, p.binxy);
end

%% Make figure
% Number of rows (first panel is for the stack-median)
nrows = ceil((nframes + 1) / p.panelsperrow);

% Median image (intial)
im_med = median(im_photon,3);

% Figure
hfig = figure('Position',[50 80 1600 900], 'Name', fign);

% Draw the median image
subplot(nrows, p.panelsperrow, 1);
hmed = imshow(im_med,[]);
title('Median')

% Draw ui control
hlist_laststate = 1 : nframes;
hlist = uicontrol('Style', 'listbox', 'String', num2cell(1:nframes), 'Max', nframes,...
    'Position', [50, 650, 70, 100], 'Value', hlist_laststate);


% Buttons
hcheck = uicontrol('Style', 'pushbutton', 'String', 'Check', 'Position',...
    [50 625 70 20], 'Callback', @FrameSelect);

hdone = uicontrol('Style', 'pushbutton', 'String', 'Done', 'Position',...
    [50 600 70 20], 'Callback', @DoneButton);

% Initialize a cell to contain individual handles
hcell = cell(nframes,1);

% Draw individual frames
for i = 1 : nframes
    % Subplot
    subplot(nrows, p.panelsperrow, i + 1);
    
    % Draw
    hcell{i} = imshow(im_photon(:,:,i),[]);
    
    title(sprintf('Frame: %i: %s', i, 'selected'));
end

%% Frame select subfunction
    function FrameSelect(~, ~)
        % Get value
        state = hlist.Value;
        
        % Which ones to remove
        removeind = ~ismember(hlist_laststate, state);
        removeind = hlist_laststate(removeind);
        
        % Which ones to add
        addind = ~ismember(state, hlist_laststate);
        addind = state(addind);
        
        % Loop through and remove
        if ~isempty(removeind)
            for j = 1 : length(removeind)
                % Get index
                fi = removeind(j);
                
                % Goto subplot and change title
                subplot(nrows, p.panelsperrow, fi + 1);
                title(sprintf('Frame: %i: %s', fi, 'unselected'));
                
                % Fade out the image
                hcell{fi}.AlphaData = 0.1;
                
            end
        end
        
        % Loop through and add
        if ~isempty(addind)
            for j = 1 : length(addind)
                % Get index
                fi = addind(j);
                
                % Goto subplot and change title
                subplot(nrows, p.panelsperrow, fi + 1);
                title(sprintf('Frame: %i: %s', fi, 'selected'));
                
                % Unfade the image
                hcell{fi}.AlphaData = 1;
            end
        end
        
        % Change the median image
        hmed.CData = median(im_photon(:,:,state),3);
        
        hlist_laststate = state;
    end

    function DoneButton(~,~)
        % Output
        F = hlist.Value;
        disp(F)        
        
        % Close figure
        close(hfig);
    end

end