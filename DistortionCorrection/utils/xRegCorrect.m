function [flimCorrected, tform, D, sbxCropped, flimShifted] = xRegCorrect(sbx, flim, getPoints)
% Inputs :
%   - sbx : sbx reference photon count image
%   - flim : photon count image obtained with FLIM dedicated sensor
%   - getPoints (default false) : boolean to determine whether to pick new
%                                 calibration points for the correction.
% Outputs :
%   - flimCorrected : corrected flim detector image
%   - tform : transformation matrix for shifting the flim images and
%             aligning them with the sbx images
%   - D : displacement field to apply to correct for FLIM sensor distortion
%   - sbxCropped : Cropped sbx image used in the function, same size as
%                  flimCorrected
%   - flimShifted : flim image aligned to sbx image prior to distortion
%                   correction
% Description :
%     This is the main function to correct the distortion of the photon
%     count image obtained with the FLIM sensor, the outputs of this
%     function can be used to then apply the distortion correction to
%     photon arrival files.
% 
if nargin < 3
    getPoints = false;
end

%% Cropping
% Cropping the images to get the same size
margin = 100;
[sbxCropped, flim] = xregCropping(sbx, flim, margin);

%% Shifting
% Aligning the flim sensor image to the sbx reference
[flimShifted, tform] = xregShift(sbxCropped, flim);

%% GUI to select points
% Allows to select new reference points to compute the distortion
% correction function. This step is only executed if the default distortion
% correction did not result in good correction.
if getPoints
    [x, y] = xregGetPoints(sbxCropped, flimShifted);
end

%% Fitting of the distortion
% Either fitting the new user defiend calibration points with the model
% function (cos^2) or loading precomputed function.
% Repeating the deformation function in a displacement field matrix.
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
% Applying the displacement field matrix to correct the distortion
flimCorrected = xregCorrectDistortion(flimShifted, D);

end
