clearvars

sbx_path = 'H:\2p\stephen\SZ333\200305_SZ333\200305_SZ333_run1\SZ333_200305_001.sbx';
flim_path = 'H:\2p\stephen\SZ333\200305_SZ333\FLIM\Output\200305_SZ333_run1_photons.tif';
tif_path = 'H:\2p\stephen\SZ333\200305_SZ333\200305_SZ333_run1\SZ333_200305_001_PMT0_binxy1_bint1_Frames50-10000.tif';
tm_path = 'H:\2p\stephen\SZ333\200305_SZ333\FLIM\Output\200305_SZ333_run1_tm.tif';

%% Loading images
cind = 1;
if ~isfile(tif_path)
    sbx2tiff_SZ(sbx_path, 0, 1, 1, [50, 10000])
end
im_sbx = readtiff(tif_path);
im_flim = readtiff(flim_path);
im_flim = im_flim(:, :, cind);

im_tm = readtiff(tm_path);
im_tm = mat2gray(im_tm(:, :, cind));

%% Taking mean image as sbx source
im_sbx = mat2gray(mean(im_sbx, 3));

%% Cross-registration
im_flim=mat2gray(im_flim);
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