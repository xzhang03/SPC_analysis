clearvars -except sbx_path flim_path
%% Loading images

disp('### Loading images')
im_sbx = Tiff(sbx_path,'r');
im_sbx = mat2gray(read(im_sbx));
im_flim = Tiff(flim_path,'r');
im_flim = mat2gray(read(im_flim));

%% Cropping
disp('### Cropping')
w = [min(size(im_sbx,1),size(im_flim,1))-100,min(size(im_sbx,2),size(im_flim,2))-100];
flim = center_patch(im_flim,w);
sbx = center_patch(im_sbx,w);

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

[optimizer,metric] = imregconfig('multimodal');
disp('## Regular registration...')
tform = imregtform(moving,fixed,rtype,optimizer, metric);
disp(tform.T)
disp('## Phase correlation...')
tformP = imregcorr(moving,fixed,rtype);
disp(tformP.T)

Rfixed = imref2d(size(fixed));
flimReg = imwarp(moving,tform,'OutputView',Rfixed);
flimRegPhase = imwarp(moving,tformP,'OutputView',Rfixed);

%% Plot output of registration
figure
subplot(311)
imshowpair(sbx,moving);
title('Initial')
subplot(312)
imshowpair(sbx,flimReg);
title('Regular registration')
subplot(313)
imshowpair(sbx,flimRegPhase);
title('Phase registration')
suptitle('Image alignment')

% %% binning pixels parce qu'on sait plus quoi faire
% sbxBinned = imresize(sbx,0.125);
% flimBinned = imresize(flimRegPhase,0.125);
% figure
% imshowpair(sbxBinned,flimBinned,'montage')
% title('Binning pixels')

%% Grid search init
w_array = linspace(600, 10000, 100);
a_array = linspace(10, 20000, 100);
[A, W] = ndgrid(a_array, w_array);
MI_grid = zeros(length(a_array), length(w_array));
CC_grid = zeros(length(a_array), length(w_array));
MSE_grid = zeros(length(a_array), length(w_array));

% Grid search
source_image = sbx;
distorted = flimRegPhase;

% cosine fitting
for IDX = 1:numel(A)
    a = A(IDX);
    w = W(IDX);
    corrected = spcDistortionCorrection(distorted, a, w);
    MI_grid(IDX) = mi(corrected(:,[1:200,end-200:end]), source_image(:,[1:200,end-200:end]));
    img_corr = corrcoef(corrected(:,[1:200,end-200:end]), source_image(:,[1:200,end-200:end]));
%     MI_grid(IDX) = mi(corrected,source_image);
%     img_corr = corrcoef(corrected,source_image);
    CC_grid(IDX) = img_corr(1, 2);
%     MSE_grid(IDX) = immse(corrected,source_image);
    MSE_grid(IDX) = immse(corrected(:,[1:200,end-200:end]), source_image(:,[1:200,end-200:end]));
end

%% Plot results
figure
subplot(121)
heatmap(round(w_array), round(a_array), MI_grid)
ylabel('Amplitude')
xlabel('Period')
title('MI')
subplot(122)
heatmap(round(w_array), round(a_array), CC_grid)
ylabel('Amplitude')
xlabel('Period')
title('CC')
% subplot(133)
% heatmap(round(w_array), round(a_array), MSE_grid)
% ylabel('Amplitude')
% xlabel('Period')
% title('MSE')

%%
flimCorrected = spcDistortionCorrection(flimRegPhase,30000,7000);

figure
% subplot(311)
% imshowpair(sbx,flim);
% title('initial')
% subplot(312)
% imshowpair(sbx,flimRegPhase);
% title('Phase registration')
% subplot(313)
imshowpair(sbx,flimCorrected);
title('Distortion corrected')

%%
plot(mean(sbx,1)-mean(flim,1))