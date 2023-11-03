function output_struct = spcSMCanalysis(smc, varargin)
% FLIM analysis based on compressed data. This method is fast and has low
% RAM requirements. SMC = sparse matrix compression

if nargin < 2
    varargin = {};
end

%% Parse inputs
p = inputParser;

% Time domain: tm and iem
addOptional(p, 'tbins', 256); % Time bins
addOptional(p, 'tcycle', 12500); % in ps

% Binning
addOptional(p, 'bint', []);

% Percentile for pctIEM (instead of using max, experimental)
addOptional(p, 'pct', 99);

% IRF
addOptional(p, 'deconvforiem', false); % If set true it takes 10x as long
addOptional(p, 'irf', [882; 6176; 15000; 15000; 3529]);

% Old parameter
addOptional(p, 'p', []);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Initialize
% Parse
timeseries = isfield(smc, 'ind2');

% Update binning
if isempty(p.bint)
    if isempty(p.p)
        % Ask
        p.bint = input('Time binning = ');
    else
        p.bint = p.p.bint;
    end
end

% Initialize
if timeseries
    % Get dimensions
    ind2 = [smc(:).ind2];
    fmax = smc(end).ind2;
else
    fmax = 1;
    ind2 = ones(size(smc));
end

% Initialize
tm = nan(fmax,1);
iem = nan(fmax,1);
photons = nan(fmax,1);
pctiem = nan(fmax,1);

%% Loop
for i = 1 : fmax

    % Current indices
    inds_curr = ind2 == i;

    % Time vectors
    if i == 1
        % Time resolution
        tres = p.tcycle / p.tbins * p.bint;
        
        % Time vector
        ivec = round(1/p.bint : sum(inds_curr))';
        tvec = (ivec - ivec(1)) * tres;

        % Initialize vtotal
        vtotal = zeros(size(ivec));
        vtotals = vtotal * ones(1,fmax);
    end
    
    % Grab temp
    smc_temp = smc(inds_curr);
    for tind = 1 : length(smc_temp)
        v = smc_temp(tind).v;
        if issparse(v)
            v = uint16(full(v))+1;
        end
        vtotal(tind) = sum(v);
    end
    vtotals(:,i) = vtotal;
    
    % Fill
    photons(i) = sum(vtotal);
    tm(i) = sum(vtotal .* tvec) / sum(vtotal);
    
    % Fill iem
    if p.deconvforiem
        vtotal = deconvlucy(vtotal, p.irf);
    end
    iem(i) = sum(vtotal) / max(vtotal) * tres;
    pctiem(i) = sum(vtotal) / prctile(vtotal, p.pct) * tres;
end

%% Output
output_struct = struct('tm', tm, 'iem', iem, 'pctiem', pctiem, 'photons', photons, 'raw', vtotals);

end