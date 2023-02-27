function [mov, params] = spcLoadsdt(fpath, rowinc, movdim)
% Load sdt file and align pixels
% [mov, params] = spcLoadsdt(fpath, rowinc, movdim)

if nargin < 3
    movdim = [511 1250];
    if nargin < 2
        rowinc = 2044.72;
    end
end

% Use bioformat toolbox to read raw data
result = bfopen(fpath);

% Get basic info
params = result{2};
l = size(result{1},1);
sz = size(result{1}{1,1});
np = sz(1) * sz(2);

%% Align
% Find indices
if isempty(movdim)
    width = floor(rowinc);
else
    width = movdim(2);
end
m = round(1 : rowinc : np - width)';
y = length(m);

% Get index array
indarray = m * ones(1,width) + ones(y,1) * (1:width) - 1;

% Initialize
if isempty(movdim)
    h = y;
else
    h = movdim(1);
end
mov = uint16(zeros(h, width, l));

%% Read in
for i = 1 : l
    f = result{1}{i,1}';
    g = f(indarray);
    
    if h == y
        mov(:,:,i) = g;
    else
        mov(:,:,i) = imresize(g, [h, width], 'bilinear');
    end
end




end