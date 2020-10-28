clearvars

% Paths for the different files
sbx_path = 'H:\2p\stephen\SZ317\SZ317_200130_001.sbx';
flim_path = 'H:\2p\stephen\SZ317\FLIM\200130_SZ317\200130_SZ317_slice1_photons.tif';
tm_path = 'H:\2p\stephen\SZ317\FLIM\200130_SZ317\200130_SZ317_slice1_tm.tif';

%% Loading images
cind = 1;
im_sbx = readtiff(sbx_path);
im_flim = readtiff(flim_path);
im_flim = im_flim(:, :, cind);
im_flim=mat2gray(im_flim);

im_tm = readtiff(tm_path);
im_tm = mat2gray(im_tm(:, :, cind));

%% Taking mean image as sbx source
im_sbx = mat2gray(mean(im_sbx, 3));

%% Cross-registration
% with getNewPoints set to false, xRegCorrect applies a preset correction
% to the photon count image acquired by the FLIM dedicated sensor.
% If this correction does not appear to match the sbx image, getNewPoints
% can be set to true, in that case the user will have to define 10
% calibration points, clicking alternatively on the 2 images on matching
% landmarks, starting with the left image.
getNewPoints = false;
[flimCorrected, tform, D, sbxCropped, flimShifted] = xRegCorrect(im_sbx, im_flim, getNewPoints);

%% Show correct FLIM photon count
figure
subplot(211)
imshowpair(sbxCropped, flimShifted);
title('Uncorrected')
subplot(212)
imshowpair(sbxCropped, flimCorrected);
title('Corrected')

%% Applying corrections to TM image
% The distortion correction can now be applied to arrival time images, the
% output will now match the sbx photon count images.
tm = crop2reference(im_tm, flimCorrected);
tmCorrected = xregCorrectTM(tm, tform, D);

%% Show corrected TM
figure
subplot(211)
imshowpair(sbxCropped, tm);
title('Uncorrected')
subplot(212)
imshowpair(sbxCropped, tmCorrected);
title('Corrected')