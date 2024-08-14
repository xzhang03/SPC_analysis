function mov = spcDecompress(smstruct, varargin)
%spcDecompress decompresses a sparse matrix structure to a lifetime decay
%matrix

%% Parse inputs
if nargin < 2
    varargin = {};
end

p = inputParser;
addOptional(p, 'type', 'uint16');

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Initialize
sz = smstruct(1).size;
sz = sz(1:4);
if sz(4) == 1
    do3d = true;
    ttind = 1;
else
    do3d = false;
end
switch p.type
    case 'uint8'
        mov = uint8(zeros(sz));
    case 'uint16'
        mov = uint16(zeros(sz));
end

%% Loop and decompress
for ind = 1 : length(smstruct)  
    % Advance index
    tind = smstruct(ind).ind;
    
    % 4d
    if ~do3d
        ttind = smstruct(ind).ind2;
    end

    % Pixels
    r = smstruct(ind).r;
    l = length(r);

    % No photon
    if l == 0
        continue;
    end

    % Load up
    c = smstruct(ind).c;
    v = smstruct(ind).v;
    if issparse(v)
        v = uint16(full(v))+1;
    end
    
    
    % Check uint8
    if strcmp(p.type, 'uint8')
        if any(v > 255) 
            fprintf('Input matrix values too large, switching to uint16\n');
        else
            v = uint8(v);
        end
    end

    for i = 1 : l
        mov(r(i),c(i),tind, ttind) = v(i);
    end
end

%% Report
% Bytes
sz = double(sz);
MB_out = sz(1) * sz(2) * sz(3) * sz(4) * 2 / 1000000;
MB_in = whos('smstruct');
MB_in = MB_in.bytes / 1000000;
fprintf('Decompressed from %0.1f MB to %0.1f MB (%0.0fx)\n', MB_in, MB_out, MB_out/MB_in);

end