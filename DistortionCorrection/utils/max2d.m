function [idx] = max2d(mat, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
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
