function [flimCorrected] = xregCorrectDistortion(sbx,flimReg, field)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
D = zeros(size(sbx,1),size(sbx,2),2);
D(:,:,1) = repmat(field,[size(sbx,1),1]);
flimCorrected = imwarp(flimReg,D,'cubic');

end

