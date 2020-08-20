function [flimCorrected,shifts,dField] = xRegCorrect(sbx,flim,getPoints)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    getPoints=false;
end

[flimReg, tform] = xregShift(sbx,flim);
shifts = tform(3,2:3);
%% GUI to select points
if getPoints
    [x,y] = xregGetPoints(sbx,flimReg);
end

%% Fitting of the distortion
[xq,approx] = xregFitDistortion(sbx,x,y);

figure
plot(x,y,'o',xq,approx,':.')
xlim([0 size(sbx,2)])

%% Correction of the distortion
flimCorrected = xregCorrectDistortion(sbx,flimReg,approx);

end

