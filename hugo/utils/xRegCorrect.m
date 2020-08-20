function [flimCorrected,shifts,D] = xRegCorrect(sbx,flim,getPoints)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    getPoints=false;
end

[flimReg, tform] = xregShift(sbx,flim);
shifts = tform.T(3,2:3);
xq = 1:size(sbx,2);
%% GUI to select points
if getPoints
    [x,y] = xregGetPoints(sbx,flimReg);
end

%% Fitting of the distortion
if getPoints
    [~,field] = xregFitDistortion(sbx,x,y);
else
    saved_model = load('modelfunction.mat');
    field = saved_model.modelfun(saved_model.beta,xq);
end
D = zeros(size(sbx,1),size(sbx,2),2);
D(:,:,1) = repmat(field,[size(sbx,1),1]);

%% Correction of the distortion
flimCorrected = xregCorrectDistortion(flimReg,D);

end

