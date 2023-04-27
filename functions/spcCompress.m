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

% 3d (single time piont) or 4d (multiple time points)
do3d = length(sz_in) == 3;
do4d = length(sz_in) == 4;

if do3d
    smstruct = struct('ind', uint16(0), 'r', uint16(0), 'c', uint16(0), 'v', uint16(0), 'size', []);
    smstruct = repmat(smstruct, [sz_in(3), 1]);
    sz_in(4) = 1;
elseif do4d
    smstruct = struct('ind2', uint16(0),'ind', uint16(0), 'r', uint16(0), 'c', uint16(0), 'v', uint16(0), 'size', []);
    smstruct = repmat(smstruct, [sz_in(3) * sz_in(4), 1]);
end

%% Loop through
ind = 0;
for ttind = 1 : sz_in(4)
    for tind = 1 : sz_in(3)
        % Advance index
        ind = ind + 1;
        
        % Fill basic
        smstruct(ind).ind = tind;
        
        if do4d
            smstruct(ind).ind2 = ttind;
        end
        
        % Get frame
        f = mov(:,:,tind, ttind);
        [r, c, v] = find(f);
        l = length(r);

        % Skip
        if l == 0
            smstruct(ind).r = [];
            smstruct(ind).c = [];
            smstruct(ind).v = [];
            smstruct(ind).size = [sz_in(1), sz_in(2), sz_in(3), sz_in(4), 0];
            continue;
        end

        % Load up non-zero pixels
        smstruct(ind).r = r;
        smstruct(ind).c = c;
        smstruct(ind).v = v;
        smstruct(ind).size = [sz_in(1), sz_in(2), sz_in(3), sz_in(4), l];
    end
end

%% Report
% Bytes
MB_in = sz_in(1) * sz_in(2) * sz_in(3) * sz_in(4) * 2 / 1000000;
MB_out = whos('smstruct');
MB_out = MB_out.bytes / 1000000;
fprintf('Compressed from %0.1f MB to %0.1f MB (-%0.1f%%)\n', MB_in, MB_out, 100-MB_out/MB_in*100);

end