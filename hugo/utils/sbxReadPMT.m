function x = sbxReadPMT(path, k, N, pmt, optolevel)

% img = sbxread(fname,k,N,pmt)
% no need to cd!
%
% Reads from frame k to k + (N - 1) in file fname
% 
% path  - the file path to .sbx file (e.g., 'xx0_000_001')
% k     - the index of the first frame to be read.  The first index is 0.
% N     - the number of consecutive frames to read starting with k.,
% optional
% pmt   - the number of the pmt, 0 for green or 1 for red, assumed to be 0
% optolevel - return a single optolevel instead of all. If passed an empty
% array, it will return all z levels of optotune

% If N>1 it returns a 4D array of size = [#pmt rows cols N] 
% If N=1 it returns a 3D array of size = [#pmt rows cols]
% If N<0 it returns an array to the end

% #pmts is the number of pmt channels being sampled (1 or 2)
% rows is the number of lines in the image
% cols is the number of pixels in each line
%
%
% The function also creates a global 'info' variable with additional
% informationi about the file

% Force a reload of the global info variables. Without this, troube arises
%clearvars -global info 
info = sbxInfo(path, true);
% Check if optotune was used, accounting for the version of scanbox being
% used
optotune_used = false;
if isfield(info, 'volscan') && info.volscan > 0, optotune_used = true; end
if ~isfield(info, 'volscan') && ~isempty(info.otwave), optotune_used = true; end

% Set to start at beginning if necessary
if nargin < 2, k = 0; end
% Set in to read the whole file if unset
if nargin < 3 || N < 0, N = info.max_idx + 1 - k; end
% Read a larger chunk if optotune was used
if nargin < 5, optolevel = []; end
if optotune_used && isinteger(optolevel), N = N*length(info.otwave); end
% Make sure that we don't search beyond the end of the file
if N > info.max_idx + 1 - k, N = info.max_idx + 1 - k; end

% Automatically set the PMT to be 0
if nargin < 4, pmt = 0; end

% Fix 0-to-1 indexing
pmt = pmt + 1;

if (isfield(info, 'fid') && info.fid ~= -1)
    try
        fseek(info.fid, k*info.nsamples, 'bof');
        x = fread(info.fid, info.nsamples/2*N, 'uint16=>uint16');
        x = reshape(x, [info.nchan info.sz(2) info.recordsPerBuffer N]);
    catch
        error('Cannot read frame. Index range likely outside of bounds.');
    end

    x = intmax('uint16') - permute(x, [1 3 2 4]);
    
    % Added by Arthur-- correct the output to a single PMT
    if info.nchan == 1
        if N > 1
            x = squeeze(x(1, :, :, :));
        else
            x = squeeze(x(1, :, :)); 
        end
    else
        if N > 1
            x = squeeze(x(pmt, :, :, :));
        else
            x = squeeze(x(pmt, :, :));
        end
    end
    
    % Check if optotune was used
    if optotune_used
        if ~isempty(optolevel)
            if nargin < 5, optolevel = 1; end
            % If optotune was used, subset the correct frames
            if isinteger(optolevel)
                optoframes = (0:length(x) + 1)*(length(info.otwave) - 1) + optolevel;
            else
                optoframes = [];
                for ol=1:length(optolevel)
                    optoframes = [optoframes ((0:length(x) + 1)*(length(info.otwave) - 1) + optolevel(ol))];
                end
            end
            optoframes = optoframes(optoframes <= size(x, 3));
            x = x(:, :, optoframes);
        end
    end
        
else
    x = [];
end