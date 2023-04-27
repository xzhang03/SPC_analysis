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
sz(3) = length(smstruct);
switch p.type
    case 'uint8'
        mov = uint8(zeros(sz));
    case 'uint16'
        mov = uint16(zeros(sz));
end

%% Loop and decompress
for tind = 1 : sz(3)
    % Pixels
    r = smstruct(tind).r;
    l = length(r);
    
    % No photon
    if l == 0
        continue;
    end
    
    % Load up
    c = smstruct(tind).c;
    v = smstruct(tind).v;
    
    % Check uint8
    if strcmp(p.type, 'uint8')
        if any(v > 255) 
            fprintf('Input matrix values too large, switching to uint16\n');
        else
            v = uint8(v);
        end
    end

    for i = 1 : l
        mov(r(i),c(i),tind) = v(i);
    end
end

end