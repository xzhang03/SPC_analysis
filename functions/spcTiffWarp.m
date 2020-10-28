function spcTiffWarp(mouse, date, varargin)
% spcTiffWarp corrects the warping artifact caused by resonant mirrors.
% This is a wrap function to usee Hugo's xRegCorrect function.

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Run indices
addOptional(p, 'sbxrun', []); % A single run number of the sbx file, used as reference (
addOptional(p, 'flimruns', []); % Run numbers of the flim files to warp

% Start with registered or unregistered data
addOptional(p, 'useregistered', true);

% Sbx parameters
addOptional(p, 'sbxframeinfo', [50 500]); % First and total number of frames to calculate reference image
addOptional(p, 'operator', @median); % What functions to use for calculating reference
addOptional(p, 'savesbxref', true); % Save sbx reference file

% Force saving
addOptional(p, 'force', false);

% Manual point sampling (shouldn't have to)
addOptional(p, 'getNewPoints', false);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% sbx reference
% Output paths
% run string
if p.sbxrun < 10
    runstr = sprintf('00%i', p.sbxrun);
elseif p.sbxrun < 100
    runstr = sprintf('0%i', p.sbxrun);
else
    runstr = num2str(p.sbxrun);
end

% Folders
datemouse = sprintf('%i_%s', date, mouse);
datemouserun = sprintf('%i_%s_run%i', date, mouse, p.sbxrun);
sbxrefpath = sprintf('\\\\%s\\data\\2p\\%s\\%s\\%s\\%s\\',...
p.server, p.user, mouse, datemouse, datemouserun);

% Filenames
sbxfn = sprintf('%s_%i_%s.sbx', mouse, date, runstr);
sbxrefname = sprintf('%i_%s_run%i_ref.tif', date, mouse, p.sbxrun);

% Check if a ref tif is already made
if exist(fullfile(sbxrefpath, sbxrefname), 'file')
    % if so, read it
    fprintf('Reading existing reference image...\n');
    reftif = readtiff(fullfile(sbxrefpath, sbxrefname));
    fprintf('Done.\n')
else
    % or read sbx
    fprintf('Generating reference image...\n');
    reftif_stack = sbxReadPMT(fullfile(sbxrefpath, sbxfn), p.sbxframeinfo(1) - 1, p.sbxframeinfo(2), 0);
    
    % Get a sample frame with median function or whatever user defines
    reftif = p.operator(reftif_stack, 3);
    
    % Write reference tif
    if p.savesbxref
        writetiff(reftif, fullfile(sbxrefpath, sbxrefname));
    end
    fprintf('Done.\n')
end


%% SPC reading and correction
% Number of FLIM runs to warp
nFLIM = length(p.flimruns);

for ii = 1 : nFLIM
    % Run index
    run = p.flimruns(ii);
    
    % Get paths
    spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
        'slice', p.slice, 'cdigit', p.cdigit);
    
    % Decide to do/redo warping for photon and tm images
    % See if warpped photon file already exist
    if p.force || ~exist(fullfile(spcpaths.fp_out, spcpaths.warptif_photons), 'file')
        dophoton = true;
    else
        % Ask if redo
        dophoton =...
                input(sprintf('%s already exists, redo? (1 = yes, 0 = no): ', spcpaths.warptif_photons)) == 1;
    end
    
    % See if warpped tm file already exist. If doing/redoing tm, will
    % always do/redo photons as well.
    if p.force || ~exist(fullfile(spcpaths.fp_out, spcpaths.warptif_tm), 'file')
        dotm = true;
        dophoton = true;
    else
        % Ask if redo
        dotm =...
                input(sprintf('%s already exists, redo? (1 = yes, 0 = no): ', spcpaths.warptif_tm)) == 1;
        dophoton = dotm;
    end
    
    % Do photons
    if dophoton
        % Read
        if p.useregistered
            photon_stack = readtiff(fullfile(spcpaths.fp_out, spcpaths.regtif_photons));
        else
            photon_stack = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_photons));
        end
        
        % N of slices
        nslices = size(photon_stack, 3);
        
        % Warp the first frame
        [flimCorrected, tform, D, ~, ~] =...
            xRegCorrect(reftif, photon_stack(:,:,1), p.getNewPoints);
        
        % Make a stack
        warptif_photons = flimCorrected;
        
        % Fill
        if nslices > 1
            warptif_photons = repmat(warptif_photons, [1 1 nslices]);
            
            for i = 2 : nslices
                cropped = crop2reference(photon_stack(:,:,i), flimCorrected);
                warptif_photons(:,:,i) = xregCorrectTM(cropped, tform, D);
            end
        end
        
        % Write
        writetiff(warptif_photons, fullfile(spcpaths.fp_out, spcpaths.warptif_photons));
    end
    
    % Do tm
    if dotm
        % Read
        if p.useregistered
            tm_stack = readtiff(fullfile(spcpaths.fp_out, spcpaths.regtif_tm));
        else
            tm_stack = readtiff(fullfile(spcpaths.fp_out, spcpaths.tif_tm));
        end
        
        % Make a stack
        warptif_tm = zeros(size(warptif_photons));
        
        % Fill
        for i = 1 : nslices
            cropped = crop2reference(tm_stack(:,:,i), flimCorrected);
            warptif_tm(:,:,i) = xregCorrectTM(cropped, tform, D);
        end
        
        % Write
        writetiff(warptif_tm, fullfile(spcpaths.fp_out, spcpaths.warptif_tm));
    end
end
end