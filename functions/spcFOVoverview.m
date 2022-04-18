function fovstruct = spcFOVoverview(mousecell, datecell, runcell, varargin)
% spcApplySbxROI generates fov-wide FLIM analysis

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'cdigit', 1); % Digits used for the "c" components in the file names (1, 2, or 3)

% Registration variables
addOptional(p, 'sourcetype', 'raw'); % Input type can be 'warped' or 'demonsreg'

% Threshold
addOptional(p, 'mintm', 1000); % Min tm
addOptional(p, 'minintensity', 2); % Min photon count

% Edge
addOptional(p, 'edge', 4); % Number of pixels to blank as edge on all sides

% Tm distribution (min max increment)
addOptional(p, 'TmDistribution', [1500 3000 1]);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Number of experiments
nexpts = length(mousecell);

% Number of trials (runs)
ntrials = 0;
for i = 1 : nexpts
    ntrials = ntrials + length(runcell{i});
end

% histogram x
histx = p.TmDistribution(1) : p.TmDistribution(3) : p.TmDistribution(2);

%% Initialize output
% Initialize
fovstruct = struct('mouse', '', 'date', '', 'run', [], 'intvec', [], 'intmean', [],...
    'tmvec', [], 'tmmean', [], 'tmdistx', [], 'tmdisty', [], 'tmdistymat', []);
fovstruct = repmat(fovstruct, [ntrials, 1]);

%% Loop through and process
itrial = 0;
hwait = waitbar(0, sprintf('Processing trials %i/%i.', itrial, ntrials));
for ii = 1 : nexpts
    % Identifiers
    mouse = mousecell{ii};
    date = datecell{ii};
    runs = runcell{ii};
    
    % Case and type
    mouse = upper(mouse);
    if ~ischar(date)
        date = num2str(date);
    end
    
    for run = runs
        % Index
        itrial = itrial + 1;
        waitbar(itrial/ntrials, hwait, sprintf('Processing trials %i/%i.', itrial, ntrials));
        
        % Get flim paths
        spcpaths = spcPath(mouse, date, run, 'server', p.server, 'user', p.user,...
            'slice', p.slice, 'cdigit', p.cdigit);

        % Image paths
        % Switch types
        switch p.sourcetype
            case 'raw'
                fp_photon = fullfile(spcpaths.fp_out, spcpaths.tif_photons);
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

        % Load photons
        usephotoninfo = false;
        if ~isempty(p.minintensity)
            if p.minintensity > 0
                usephotoninfo = true;
                [im_photon, ~] = readtiff(fp_photon);
                
                if p.edge > 0
                    im_photon = im_photon(p.edge+1 : end-p.edge, p.edge+1 : end-p.edge, :);
                end
            end
        end
        
        % Load Tm
        [im_tm, ~] = readtiff(fp_tm);  
        if p.edge > 0
            im_tm = im_tm(p.edge+1 : end-p.edge, p.edge+1 : end-p.edge, :);
        end
        
        % Number of reads
        nreads = size(im_tm, 3);
        
        % Initialize vecs
        intvec = zeros(nreads, 1);
        tmvec = zeros(nreads, 1);
        histy = zeros(length(histx), nreads);
        
        for iread = 1 : nreads
            if usephotoninfo
                % Mask
                fpts = im_photon(:,:,iread);
                mask = fpts >= p.minintensity;
                tms = im_tm(:,:,iread) .* mask;
                
                % Get photon info
                pts = fpts(:); % No threshold
%                 pts = fpts(mask); % If thresholded
                intvec(iread) = mean(pts);
                
            else
                tms = im_tm(:,:,iread);
            end
            
            % Threshold by minimal lifetime
            tms = tms(tms >= p.mintm);
            tmvec(iread) = mean(tms);
            histy(:,iread) = hist(tms, histx);
        end
        
        % Put info in
        fovstruct(itrial).mouse = mouse;
        fovstruct(itrial).date = date;
        fovstruct(itrial).run = run;
        fovstruct(itrial).intvec = intvec;
        fovstruct(itrial).intmean = nanmean(intvec);
        fovstruct(itrial).tmvec = tmvec;
        fovstruct(itrial).tmmean = nanmean(tmvec);
        fovstruct(itrial).tmdistx = histx;
        fovstruct(itrial).tmdisty = mean(histy,2);
        fovstruct(itrial).tmdistymat = histy;
    end
end
close(hwait)


end