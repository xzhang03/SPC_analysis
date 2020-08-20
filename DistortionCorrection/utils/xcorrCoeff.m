function [CC] = xcorrCoeff(im1, im2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
img_corr = corrcoef(im1, im2);
CC = img_corr(1, 2);
end
