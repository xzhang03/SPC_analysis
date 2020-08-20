clearvars -except sbx_path flim_path x y

sbx_path = 'H:\2p\stephen\SZ317\SZ317_200130_001.sbx';
flim_path = 'H:\2p\stephen\SZ317\FLIM\200130_SZ317\200130_SZ317_slice1_photons_reg.tif';
tif_path = 'H:\2p\stephen\SZ317\SZ317_200130_001_PMT0_binxy1_bint1_Frames50-5000.tif';
%% Loading images
if ~isfile(tif_path)
    sbx2tiff_SZ(sbx_path, 0, 1, 1, [50 5000])
end
im_sbx = readtiff(tif_path);
im_flim = Tiff(flim_path,'r'); im_flim = mat2gray(read(im_flim));

%%
im_sbx=mat2gray(mean(im_sbx,3));
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