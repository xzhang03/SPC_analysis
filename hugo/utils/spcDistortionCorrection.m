function [corrected_tm,D] = spcDistortionCorrection(unregistered_tm,a,w)
% Performs sine unwarping to correct for laser-scanning induced distortion

maxI = max(unregistered_tm(:));

yvect = 1:size(unregistered_tm,2);
proj = a * sin(1/w * (yvect - size(unregistered_tm,2)/2));

D = zeros(size(unregistered_tm,1),size(unregistered_tm,2),2);
D(:,:,1) = repmat(proj,[size(unregistered_tm,1),1]);

corrected_tm = imwarp(unregistered_tm,D,'cubic');
corrected_tm(corrected_tm>maxI) = maxI;
end

