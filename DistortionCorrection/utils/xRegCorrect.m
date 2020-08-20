function [flimCorrected, tform, D, sbxCropped, flimShifted] = xRegCorrect(sbx, flim, getPoints)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    getPoints = false;
end

%% Cropping
margin = 100;
[sbxCropped, flim] = xregCropping(sbx, flim, margin);

%% Shifting
[flimShifted, tform] = xregShift(sbxCropped, flim);
shifts = tform.T(3, 2:3);

%% GUI to select points
if getPoints
    [x, y] = xregGetPoints(sbxCropped, flimShifted);
end

%% Fitting of the distortion
xq = 1:size(sbxCropped, 2);

if getPoints
    [~, field] = xregFitDistortion(sbxCropped, x, y);
else
    saved_model = load('modelfunction.mat');
    field = saved_model.modelfun(saved_model.beta, xq);
end
D = zeros(size(sbxCropped, 1), size(sbxCropped, 2), 2);
D(:, :, 1) = repmat(field, [size(sbxCropped, 1), 1]);

%% Correction of the distortion
flimCorrected = xregCorrectDistortion(flimShifted, D);

end
