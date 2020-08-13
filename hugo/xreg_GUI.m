clearvars -except sbx_path flim_path x y

%% Loading images
im_sbx = Tiff(sbx_path,'r'); im_sbx = mat2gray(read(im_sbx));
im_flim = Tiff(flim_path,'r'); im_flim = mat2gray(read(im_flim));

%% Cropping
[sbx,flim] = xregCropping(im_sbx,im_flim);

%% Registration - Translation
flimReg = xregShift(sbx,flim);

%% GUI to select points
if ~exist('y','var')
    [x,y] = xregGetPoints(sbx,flimReg);
end

%% Fitting of the distortion
[xq,approx] = xregFitDistortion(sbx,x,y);

figure
plot(x,y,'o',xq,approx,':.')
xlim([0 size(sbx,2)])

%% Correction of the distortion
flimCorrected = xregCorrectDistortion(sbx,flimReg,approx);

figure
imshowpair(sbx,flimCorrected);
title('Distortion corrected')
