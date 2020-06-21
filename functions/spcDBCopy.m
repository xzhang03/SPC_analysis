function spcDBCopy(mouse, date, varargin)
% spcDBCopy copies key files to a dropbox folder

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Destination path
addOptional(p, 'destfp', '');

% Start with registered or unregistered data
addOptional(p, 'useregistered', true);

% Which sections to use (if not specified, it will be corrected to 2 : end
% below. If the section is deemed blurry, it will be thrown out as well
addOptional(p, 'sections', []);

% Copy local-normalized images
addOptional(p, 'copylnimages', false);

% Multiple fovs (affects cross run names)
addOptional(p, 'multifov', false);
addOptional(p, 'fov', 1); % Specify fov

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

% Run or Slice
if p.slice
    RunOrSlice = 'slice';
else
    RunOrSlice = 'run';
end

% Halt and catch fire if no destination path
if isempty(p.destfp)
    disp('Need a destination file path.')
    return;
else
    destfp = p.destfp;
end

%% IO
% Get paths (run number does not matter)
spcpaths = spcPath(mouse, date, 0, 'server', p.server, 'user', p.user,...
    'slice', p.slice, 'cdigit', p.cdigit, 'multifov', p.multifov, 'fov', p.fov);

% Load
load(fullfile(spcpaths.fp_out, spcpaths.xrun_mat), 'ROI_struct', 'nsections');

% Throw out the first section if
% p.sections is unspecified (First section contains frames where the
% resonant mirrows have yet to be up to speed)
if isempty(p.sections)
    p.sections = 2 : nsections;
end

% Get runs
runs = [ROI_struct.run];

%% Destination path
% Make mouse folder if there is none
if exist(destfp, 'dir') && ~exist(fullfile(destfp,mouse), 'dir')
    mkdir(fullfile(destfp,mouse));
end

% Make day folder if there is none
datefolder = sprintf('%s_%s', date, mouse);
destfp_full = fullfile(destfp,mouse,datefolder);
if exist(fullfile(destfp,mouse), 'dir') && ~exist(destfp_full, 'dir')
    mkdir(destfp_full);
end

%% Copy files
% Copy REF ROI image
copyfile(fullfile(spcpaths.fp_out, spcpaths.ROI_ref),...
    fullfile(destfp_full, spcpaths.ROI_ref))

% Copy REF neuropil ROI image
copyfile(fullfile(spcpaths.fp_out, spcpaths.ROInp_ref),...
    fullfile(destfp_full, spcpaths.ROInp_ref))

% Copy tm csv file
copyfile(fullfile(spcpaths.fp_out, spcpaths.xruntm_csv),...
    fullfile(destfp_full, spcpaths.xruntm_csv))

% Copy tm csv file
copyfile(fullfile(spcpaths.fp_out, spcpaths.xruntmnp_csv),...
    fullfile(destfp_full, spcpaths.xruntmnp_csv))

% Copy crossrun mat file
copyfile(fullfile(spcpaths.fp_out, spcpaths.xrun_mat),...
    fullfile(destfp_full, spcpaths.xrun_mat))

%% Make image files
% Loop through runs
for ind = 1 : length(runs)
    % Get run
    run = runs(ind);
    
    % Image filename (output)
    fn_tm_out = sprintf('%s_%s_%s%i_tm_med.tif', date, mouse, RunOrSlice, run);
    fn_photons_out = sprintf('%s_%s_%s%i_photons_med.tif', date, mouse, RunOrSlice, ...
    run);

    % Image filename (input)
    if p.useregistered
        fn_tm_in = sprintf('%s_%s_%s%i_tm_reg.tif', date, mouse, RunOrSlice, run);
        fn_photons_in = sprintf('%s_%s_%s%i_photons_reg.tif', date, mouse, RunOrSlice, ...
    run);
    else
        fn_tm_in = sprintf('%s_%s_%s%i_tm.tif', date, mouse, RunOrSlice, run);
        fn_photons_in = sprintf('%s_%s_%s%i_photons.tif', date, mouse, RunOrSlice, ...
    run);
    end
    
    % Sections
    % Find useful sections (not blurry or specified as not used here)
    sections_to_use = intersect(p.sections, ROI_struct(ind).sections);
    
    % Read images
    tm_stack = readtiff(fullfile(spcpaths.fp_out, fn_tm_in));
    photons_stack = readtiff(fullfile(spcpaths.fp_out, fn_photons_in));
    
    % Take median
    tm_stack_med = median(tm_stack(:,:,sections_to_use), 3);
    photons_stack_med = median(photons_stack(:,:,sections_to_use), 3);
    
    % Write tiffs
    writetiff(tm_stack_med, fullfile(destfp_full,fn_tm_out), 'double');
    writetiff(photons_stack_med, fullfile(destfp_full,fn_photons_out), 'double');
    
    % Write local-normalized images
    if p.copylnimages
        fn_ln = sprintf('%s_%s_%s%i_photons_ln.tif', date, mouse, RunOrSlice, ...
    run);
        writetiff(imresize(ROI_struct(ind).im,4), fullfile(destfp_full,fn_ln), 'double');
    end
   
end
end