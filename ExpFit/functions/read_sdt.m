function [phot_count,out_image] = read_sdt(path)%,crop)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% if nargin < 2
%     crop = true;
% end

sdt = bfopen(path);
out_image = cat(3,sdt{1}{:,1});

% if crop
%     temp = zeros(size(out_image),'like',out_image);
%     for layer = 1:size(out_image,3)
%         temp(:,:,layer) = 
phot_count = squeeze(sum(sum(out_image,1),2));

end

