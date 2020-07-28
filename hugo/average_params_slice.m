clear

folder = 'H:\2p\stephen\SZ317\FLIM\200130_SZ317\200130_SZ317_slice1';
fn_photons_in = '200130_SZ317_slice1_c%02d_photons.asc';
fn_tm_in = '200130_SZ317_slice1_c%02d_t1.asc';
cinds = 1:90;

%% Loading photon and tm files
w_array = logspace(2,5,50);
a_array = logspace(2,5,50);
[A,W] = ndgrid(a_array,w_array);
MI_grid_all = zeros(length(w_array),length(a_array),numel(cinds));
corrCoeff_grid_all = zeros(length(w_array),length(a_array),numel(cinds));

%%

parfor ID = cinds
    im_photon = load(fullfile(folder, sprintf(fn_photons_in,ID)));

    im_tm = load(fullfile(folder, sprintf(fn_tm_in,ID)));

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
    for IDX = 1:numel(A)
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

%%
figure
subplot(121)
imshow(mean(MI_grid_all(:,:,1),3))
colormap(gca,'parula')
colorbar
title('Mutual information')
subplot(122)
imshow(mean(corrCoeff_grid_all(:,:,1),3))
colorbar

colormap(gca,'parula')
title('Correlation coefficient')

%%
figure
for i = 1:10
    subplot(2,5,i)
    imshow(corrCoeff_grid_all(:,:,i))
    colormap(gca,'parula')
    colorbar
end
%% Check idx of best MI for all images

best_fit_MI = zeros(numel(cinds),1);
best_fit_corr = zeros(numel(cinds),1);
for ID = cinds
    loc = max2d(MI_grid_all(:,:,ID));
    best_fit_MI(ID) = loc;
    loc = max2d(abs(corrCoeff_grid_all(:,:,ID)));
    best_fit_corr(ID) = loc;
end


