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

[xq,approx,beta,modelfun] = xregFitDistortion(sbx,x,y);


figure
plot(x,y,'o',xq,approx,':.')
xlim([0 size(sbx,2)])

%% Correction of the distortion
flimCorrected = xregCorrectDistortion(sbx,flimReg,approx);

figure
imshowpair(sbx,flimCorrected);
title('Distortion corrected')

%%
figure
subplot(211)
imshowpair(sbx,flim,'montage')
title('SBX PMT photon image and FLIM PMT photon image')
subplot(212)
imshowpair(sbx,flimReg)
title('Overlayed')

%%
figure
plot(x,y,'o',xq,approx,':.')
xlabel('Horizontal position [pixels]')
ylabel('Distortion amplitude [pixels]')
xlim([0 size(sbx,2)])
title('Fit of the measured distortion')
%%
figure
imshowpair(sbx,flimCorrected)
title('Corrected image overlayed on SBX PMT image')
