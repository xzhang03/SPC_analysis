function frames = spcSMCphotonframe(smc, varargin)

if nargin < 2
    varargin = {};
end

%% Parse inputs
p = inputParser;

% Frame range
addOptional(p, 'frames', []); % Leave empty to use all
addOptional(p, 'combined', false); % Combine all the frames into 1 if true

% Crop
addOptional(p, 'crop', []);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Initialize
% Parse
timeseries = isfield(smc, 'ind2');

% Get size
sz = smc(1).size;
sz = sz(1:4);

% Calculate the number of frames
if p.combined
    nframes = 1;
elseif ~isempty(p.frames)
    nframes = length(p.frames);
else
    nframes = sz(4);
end

% Initialize mov
frames = uint16(zeros(sz(1), sz(2), nframes));

% Initialize
if timeseries
    % Get dimensions
    ind2 = [smc(:).ind2];
    fmax = smc(end).ind2;
else
    fmax = 1;
    ind2 = ones(size(smc));
end


%% Loop
ind = 0;
for i = 1 : fmax
    % Check if we care about this frame
    if ~isempty(p.frames)
        if ~ismember(i, p.frames)
            % Doesn't belong to any of the designated frames, skip
            continue;
        end
    end
    
    % Get frame
    if p.combined
        ind = 1;
    else
        ind = ind + 1;
    end
    f = frames(:,:,ind);

    % Write pixels
    inds_curr = ind2 == i;
    smc_temp = smc(inds_curr);
    for tind = 1 : length(smc_temp)
        % Get vectors
        r = smc_temp(tind).r;
        c = smc_temp(tind).c;
        v = smc_temp(tind).v;
        if issparse(v)
            v = uint16(full(v))+1;
        end

        for pind = 1 : length(r)
            f(r(pind), c(pind)) =  f(r(pind), c(pind)) + uint16(v(pind));
        end
    end

    % Put frames back
    frames(:,:,ind) = f;
end

%% Crop
if ~isempty(p.crop)
    frames = frames(p.crop(2):p.crop(4), p.crop(1):p.crop(3), :);
end
end
