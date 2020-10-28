function binned = bin2d(stack, factor)
%BINXY Adapted from Rohan's binning script, bins only in XY

    if nargin < 2, factor = 2; end
    factor = round(factor);

    [y, x, nframes] = size(stack);
    
    % Remove edge pixels if necessary
    if mod(x, factor), x = x - mod(x, factor); end
    if mod(y, factor), y = y - mod(y, factor); end
    stack = stack(1:y, 1:x, :);

    % Turn movie into 2-D vector
    stack = reshape(stack, y, x*nframes);
    [m, n] = size(stack);

    % Bin along columns:
    stack = sum(reshape(stack, factor, []), 1);

    % Bin along rows:
    stack = reshape(stack, m/factor, []).'; %Note transpose
    stack = sum(reshape(stack, factor, []), 1);
    stack = reshape(stack, n/factor, []).'; %Note transpose 

    % Turn back into original shape:
    binned = reshape(stack, y/factor, x/factor, nframes);
end

