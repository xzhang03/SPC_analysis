function [idx] = max2d(mat, varargin)
% Returns the coordinates of the max value of a 2d array.
if nargin < 2
    form = '2d';
else
    form = varargin(1);
end

[~, idx] = max(mat(:));

if strcmp(form, '2d')
    [I, J] = ind2sub(size(mat), idx);
    idx = [I, J];
end
end
