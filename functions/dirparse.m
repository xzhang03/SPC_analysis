function outputpath_struct = dirparse(pathin, nlevels)
% fpathparse parses path to get different levels of parent paths (windows)
% outputpath_struct = dirparse(pathin, levels)

if nargin < 2
    nlevels  = 3;
end

% Check if the input path ends with a \
if strcmpi(pathin(end), '\')
    pathin = pathin(1:end-1);
end

% Find the levels
inds = strfind(pathin, '\');

if length(inds) >= nlevels
    % More levels than asked
    inds = inds(end:-1:end-nlevels+1);
else
    % Fewer levesl than asked
    nlevels = length(inds);
    inds = inds(end:-1:1);
end

% Output
outputpath_struct = struct('path', '');
outputpath_struct = repmat(outputpath_struct, [nlevels 1]);

for i = 1 : nlevels
    outputpath_struct(i).path = pathin(1:inds(i)-1);
end

end