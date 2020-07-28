clear
mouse = 'SZ334';
date = 200624;
run = 1;
user = 'stephen';

%% Get paths
spcpaths = spcPath(mouse, date, run, 'user',user);
inds = spcpaths.cinds;

%% Loading photon and tm files
w_array = 300:1000:200000;
a_array = 1:500:100000;
[A,W] = ndgrid(a_array,w_array);
MI_grid_all = zeros(length(w_array),length(a_array),numel(inds));
corrCoeff_grid_all = zeros(length(w_array),length(a_array),numel(inds));

%%
for ID = inds
    im_photon = load(fullfile(spcpaths.fp, sprintf(spcpaths.photons_in,ID)));

    im_tm = load(fullfile(spcpaths.fp, sprintf(spcpaths.tm_in,ID)));

    tm = mat2gray(im_tm);
    phot = mat2gray(im_photon);

    % Mutual information before any operations
    MI_initial = mi(tm,phot);

    % Normalized x-correlation to get shifts
    xC_initial = normxcorr2(tm,phot);
    [~, ind] = max(xC_initial(:));
    [X, Y] = ind2sub(size(xC_initial),ind);
    shifts = [size(tm,1)-X, size(tm,2)-Y];
    % figure
    % surf(xC_initial)
    % shading flat
    MI_grid = zeros(length(w_array),length(a_array));
    
    corrCoeff_grid = zeros(length(w_array),length(a_array));
    % cosine fitting
    parfor IDX = 1:numel(A)
        a = A(IDX);
        w = W(IDX);
        corrected = spcDistortionCorrection(tm,a,w);
        MI_grid(IDX) = mi(corrected, phot);
        img_corr = corrcoef(corrected,phot);
        corrCoeff_grid(IDX) = img_corr(1,2);
    end
    MI_grid_all(:,:,ID) = MI_grid;
    corrCoeff_grid_all(:,:,ID) = corrCoeff_grid;
end

% %
% outgrid = mat2gray(MI_grid);
% outgrid = imgaussfilt(outgrid,20);
% figure
% imshow(outgrid')

%% Check idx of best MI for all images

best_fit_MI = zeros(numel(inds),1);
best_fit_corr = zeros(numel(inds),1);
for ID = inds
    loc = max2d(MI_grid_all(:,:,ID));
    best_fit_MI(ID) = loc;
    loc = max2d(corrCoeff_grid_all(:,:,ID));
    best_fit_corr(ID) = loc;
end

%%
a = 200;
w = 3500;

im_photon = load(fullfile(spcpaths.fp, sprintf(spcpaths.photons_in,1)));

im_tm = load(fullfile(spcpaths.fp, sprintf(spcpaths.tm_in,1)));

tm = mat2gray(im_tm);
phot = mat2gray(im_photon);
corrected = spcDistortionCorrection(tm,a,w);
MI_initial = mi(tm,phot);
CC_initial = corrcoef(tm,phot);
CC_initial = CC_initial(1,2);
MI = mi(corrected,phot);
CC = corrcoef(corrected,phot);
CC = CC(1,2);

figure
subplot(311)
imshow(tm)
title(['Initial TM, ', 'MI : ',num2str(MI_initial), ', CC : ',num2str(CC_initial)])
subplot(312)
imshow(phot)
title('Photon')
subplot(313)
imshow(corrected)
title(['Corrected TM, ', 'MI : ',num2str(MI), ', CC : ',num2str(CC)])
