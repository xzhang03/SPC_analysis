clearvars -except sbx_path flim_path x y

%% Loading images
im_sbx = Tiff(sbx_path,'r'); im_sbx = mat2gray(read(im_sbx));
im_flim = Tiff(flim_path,'r'); im_flim = mat2gray(read(im_flim));

%% Cropping
[sbx,flim] = xregCropping(im_sbx,im_flim);

%%
[flimCorrected,shifts,D] = xRegCorrect(sbx,flim,false);

%%
figure
imshowpair(sbx,flimCorrected);
title('Distortion corrected')
