function smstruct = spcCompress(mov, varargin)
%spcCompress compresses sdt file into a sparse matrix structure

%% Parse inputs
if nargin < 2
    varargin = {};
end

p = inputParser;
addOptional(p, 'temp', true);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Initialize
sz_in = size(mov);
smstruct = struct('ind', uint8(0), 'r', uint16(0), 'c', uint16(0), 'v', uint16(0), 'size', []);
smstruct = repmat(smstruct, [sz_in(3), 1]);

%% Loop through
for tind = 1 : sz_in(3)
    % Fill basic
    smstruct(tind).ind = tind;

    % Get frame
    f = mov(:,:,tind);
    [r, c, v] = find(f);
    l = length(r);
    
    % Skip
    if l == 0
        smstruct(tind).r = [];
        smstruct(tind).c = [];
        smstruct(tind).v = [];
        smstruct(tind).size = [sz_in(1), sz_in(2), 0];
        continue;
    end

    % Load up non-zero pixels
    smstruct(tind).r = r;
    smstruct(tind).c = c;
    smstruct(tind).v = v;
    smstruct(tind).size = [sz_in(1), sz_in(2), l];
end

%% Report
% Bytes
MB_in = sz_in(1) * sz_in(2) * sz_in(3) * 2 / 1000000;
MB_out = whos('smstruct');
MB_out = MB_out.bytes / 1000000;
fprintf('Compressed from %0.1f MB to %0.1f MB (-%0.1f%%)\n', MB_in, MB_out, 100-MB_out/MB_in*100);

end