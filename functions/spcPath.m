function spcpaths = spcPath(mouse, date, run, varargin)
% spcPath generates all the path info for FLIM data
%% Parse inputs
p = inputParser;

addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

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

%% Input files
% Run or Slice
if p.slice
    RunOrSlice = 'slice';
else
    RunOrSlice = 'run';
end

% File path
spcpaths.fp = sprintf('\\\\%s\\data\\2p\\%s\\%s\\%s_%s\\FLIM\\%s_%s_%s%i', p.server, ...
    p.user, mouse, date, mouse, date, mouse, RunOrSlice, run);

% Names based on the number of c digits
cstring = sprintf('c%%0%id', p.cdigit);

% Tm input file
spcpaths.tm_in = sprintf('%s_%s_%s%i_%s_t1.asc', date, mouse, RunOrSlice, run,...
    cstring);

% Photons input file
spcpaths.photons_in = sprintf('%s_%s_%s%i_%s_photons.asc', date, mouse, RunOrSlice, run,...
    cstring);

%% Get the c indices
% Dir
flist = dir(fullfile(spcpaths.fp,'*.img'));

% Get the indices
spcpaths.n = size(flist,1);
spcpaths.cinds = 1 : spcpaths.n;

%% Output files
% File path
fp_parent = sprintf('\\\\%s\\data\\2p\\%s\\%s\\%s_%s\\FLIM', p.server, ...
    p.user, mouse, date, mouse);
spcpaths.fp_out = fullfile(fp_parent, 'Output');

% Make output folder if it doesn't exist
if exist(fp_parent, 'dir') && ~exist(spcpaths.fp_out, 'dir')
    mkdir(spcpaths.fp_out);
end

% Tm output files
spcpaths.tif_tm = sprintf('%s_%s_%s%i_tm.tif', date, mouse, RunOrSlice, ...
    run);
spcpaths.regtif_tm = sprintf('%s_%s_%s%i_tm_reg.tif', date, mouse, RunOrSlice, ...
    run);

% Photon output files
spcpaths.tif_photons = sprintf('%s_%s_%s%i_photons.tif', date, mouse, RunOrSlice, ...
    run);
spcpaths.regtif_photons = sprintf('%s_%s_%s%i_photons_reg.tif', date, mouse, RunOrSlice, ...
    run);

% Mat output files
spcpaths.mat = sprintf('%s_%s_%s%i_output.mat', date, mouse, RunOrSlice, ...
    run);

%% Cross-expt output file
if p.multifov
    % Multiple fovs
    % Mat cross expt output file
    spcpaths.xrun_mat = sprintf('%s_%s_crossrun_output_fov%i.mat', date, mouse, p.fov);

    % Ref image output file
    spcpaths.ROI_ref = sprintf('%s_%s_ROIRefs_fov%i.tif', date, mouse, p.fov);

    % Csv tm output file
    spcpaths.xruntm_csv = sprintf('%s_%s_tm_fov%i.csv', date, mouse, p.fov);
    
    % Ref image neuropil output file
    spcpaths.ROInp_ref = sprintf('%s_%s_npROIRefs_fov%i.tif', date, mouse, p.fov);
    
    % Csv tm neuropil output file
    spcpaths.xruntmnp_csv = sprintf('%s_%s_nptm_fov%i.csv', date, mouse, p.fov);
else
    % Single fov
    % Mat cross expt output file
    spcpaths.xrun_mat = sprintf('%s_%s_crossrun_output.mat', date, mouse);

    % Ref image output file
    spcpaths.ROI_ref = sprintf('%s_%s_ROIRefs.tif', date, mouse);

    % Csv tm output file
    spcpaths.xruntm_csv = sprintf('%s_%s_tm.csv', date, mouse);
    
    % Ref image neuropil output file
    spcpaths.ROInp_ref = sprintf('%s_%s_npROIRefs.tif', date, mouse);

    % Csv tm neuropil output file
    spcpaths.xruntmnp_csv = sprintf('%s_%s_nptm.csv', date, mouse);
end

end