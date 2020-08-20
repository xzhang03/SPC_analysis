function [corrected_tm, D] = spcDistortionCorrection(unregistered_tm, a, w, mode)
% Performs sine unwarping to correct for laser-scanning induced distortion
if nargin < 4
    mode = 'cos';
end


yvect = 1:size(unregistered_tm, 2);

if strcmp(mode, 'cos')
    proj = a * (cos(1 / w * (yvect - size(unregistered_tm, 2) / 2)) - 1);
    proj(ceil(numel(yvect) / 2):end) = -proj(ceil(numel(yvect) / 2):end);
elseif strcmp(mode, 'sin')
    proj = a * sin(1/w*(yvect - size(unregistered_tm, 2) / 2));
end

D = zeros(size(unregistered_tm, 1), size(unregistered_tm, 2), 2);
D(:, :, 1) = repmat(proj, [size(unregistered_tm, 1), 1]);

corrected_tm = imwarp(unregistered_tm, D, 'cubic');

% Make sure that no pixel values are greater than unregistered image max.
maxI = max(unregistered_tm(:));
corrected_tm(corrected_tm > maxI) = maxI;
end
