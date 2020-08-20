clearvars -except sbx_path flim_path

%% Loading images

disp('### Loading images')
im_sbx = Tiff(sbx_path, 'r');
im_sbx = mat2gray(read(im_sbx));
im_flim = Tiff(flim_path, 'r');
im_flim = mat2gray(read(im_flim));

%% Cropping
disp('### Cropping')
w = [min(size(im_sbx, 1), size(im_flim, 1)) - 100, min(size(im_sbx, 2), size(im_flim, 2)) - 100];
flim = center_patch(im_flim, w);
sbx = center_patch(im_sbx, w);

% %% Histogram matching
% disp('### Histogram matching')
% figure
% subplot(211)
% imshowpair(sbx,flim,'montage')
% title('Initial')
%
% flimMatched = imhistmatch(flim,sbx);
% subplot(212)
% imshowpair(sbx,flimMatched,'montage')
% title('Histogram matched')

%% Registration
rtype = 'translation';
moving = flim;
fixed = sbx;

[optimizer, metric] = imregconfig('multimodal');
disp('## Regular registration...')
tform = imregtform(moving, fixed, rtype, optimizer, metric);
disp('    Shifts : ')
disp(tform.T(3, 1:2))
disp('## Phase correlation...')
tformP = imregcorr(moving, fixed, rtype);
disp('    Shifts : ')
disp(tformP.T(3, 1:2))

Rfixed = imref2d(size(fixed));
flimReg = imwarp(moving, tform, 'OutputView', Rfixed);
flimRegPhase = imwarp(moving, tformP, 'OutputView', Rfixed);

%% Plot output of registration
figure
subplot(311)
imshowpair(sbx, moving);
title('Initial')
subplot(312)
imshowpair(sbx, flimReg);
title('Regular registration')
subplot(313)
imshowpair(sbx, flimRegPhase);
title('Phase registration')

%% binning pixels parce qu'on sait plus quoi faire
sbxBinned = imresize(sbx, 0.125);
flimBinned = imresize(flimRegPhase, 0.125);
sbxBinned = imresize(sbxBinned, 8);
flimBinned = imresize(flimBinned, 8);
figure
imshowpair(sbxBinned, flimBinned, 'montage')
title('Binning pixels')

%% Grid search init
w_array = linspace(600, 10000, 50);
a_array = linspace(10, 3000, 50);
[A, W] = ndgrid(a_array, w_array);
MI_grid = zeros(length(a_array), length(w_array));
CC_grid = zeros(length(a_array), length(w_array));
MSE_grid = zeros(length(a_array), length(w_array));

MI_grid_sin = zeros(length(a_array), length(w_array));
CC_grid_sin = zeros(length(a_array), length(w_array));
MSE_grid_sin = zeros(length(a_array), length(w_array));

%% Grid search
source_image = sbx;
distorted = flimRegPhase;

% cosine fitting
parfor IDX = 1:numel(A)
    a = A(IDX);
    w = W(IDX);
    corrected = spcDistortionCorrection(distorted, a, w);
    MI_grid(IDX) = mi(corrected, source_image);
    img_corr = corrcoef(corrected, source_image);
    CC_grid(IDX) = img_corr(1, 2);
    MSE_grid(IDX) = immse(corrected, source_image);

    corrected = spcDistortionCorrection(distorted, a, w, 'sin');
    MI_grid_sin(IDX) = mi(corrected, source_image);
    img_corr = corrcoef(corrected, source_image);
    CC_grid_sin(IDX) = img_corr(1, 2);
    MSE_grid_sin(IDX) = immse(corrected, source_image);
end

%% Plot results
figure
subplot(121)
heatmap(round(w_array), round(a_array), MI_grid_sin);
ylabel('Amplitude')
xlabel('Period')
title('MI')
subplot(122)
heatmap(round(w_array), round(a_array), CC_grid_sin);
ylabel('Amplitude')
xlabel('Period')
title('CC')
% subplot(133)
% heatmap(round(w_array), round(a_array), MSE_grid)
% ylabel('Amplitude')
% xlabel('Period')
% title('MSE')

%%
flimCorrected = spcDistortionCorrection(flimRegPhase, 3000, 800);

figure
% subplot(311)
% imshowpair(sbx,flim);
% title('initial')
% subplot(312)
% imshowpair(sbx,flimRegPhase);
% title('Phase registration')
% subplot(313)
imshowpair(sbx, flimCorrected);
title('Distortion corrected')

%%
plot(mean(sbx, 1)-mean(flim, 1))