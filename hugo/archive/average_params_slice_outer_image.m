clear

folder = 'H:\2p\stephen\SZ317\FLIM\200130_SZ317\200130_SZ317_slice1';
fn_photons_in = '200130_SZ317_slice1_c%02d_photons.asc';
fn_tm_in = '200130_SZ317_slice1_c%02d_t1.asc';
cinds = 1:90;

%% Loading photon and tm files
w_array = linspace(600, 5000, 25);
a_array = linspace(1, 1500, 25);
[A, W] = ndgrid(a_array, w_array);
MI_grid_all = zeros(length(w_array), length(a_array), numel(cinds));
corrCoeff_grid_all = zeros(length(w_array), length(a_array), numel(cinds));

%

for ID = cinds
    im_photon = load(fullfile(folder, sprintf(fn_photons_in, ID)));

    im_tm = load(fullfile(folder, sprintf(fn_tm_in, ID)));

    tm = mat2gray(im_tm);
    phot = mat2gray(im_photon);

    % Mutual information before any operations
    MI_initial = mi(tm, phot);

    MI_grid = zeros(length(w_array), length(a_array));
    corrCoeff_grid = zeros(length(w_array), length(a_array));

    % cosine fitting
    parfor IDX = 1:numel(A)
        a = A(IDX);
        w = W(IDX);
        corrected = spcDistortionCorrection(tm, a, w);
        %         third = round(size(phot,2)/3);
        %         corrected_cropped = corrected(:,[1:third, 2*third:size(tm,2)]);
        %         phot_cropped =phot(:,[1:third, 2*third:size(tm,2)]);
        MI_grid(IDX) = mi(corrected, phot);
        img_corr = corrcoef(corrected, phot);
        corrCoeff_grid(IDX) = img_corr(1, 2);
    end
    MI_grid_all(:, :, ID) = MI_grid;
    corrCoeff_grid_all(:, :, ID) = corrCoeff_grid;
end

% %
% outgrid = mat2gray(MI_grid);
% outgrid = imgaussfilt(outgrid,20);
% figure
% imshow(outgrid')

%% put data in table for heatmap
MI_mean = mean(MI_grid_all, 3);
CC_mean = mean(corrCoeff_grid_all, 3);

figure
subplot(121)
heatmap(round(a_array), round(w_array), MI_mean);
subplot(122)
heatmap(round(a_array), round(w_array), CC_mean);

%%
surf(mean(MI_grid_all, 3))

%%
figure
subplot(121)
imshow(mean(MI_grid_all(:, :, 1), 3))
colormap(gca, 'parula')
colorbar
title('Mutual information')
subplot(122)
imshow(mean(corrCoeff_grid_all(:, :, 1), 3))
colorbar

colormap(gca, 'parula')
title('Correlation coefficient')

%%
figure
for i = 1:2
    subplot(2, 1, i)
    imshow(corrCoeff_grid_all(:, :, i))
    colormap(gca, 'parula')
    colorbar
end

%% Check idx of best MI for all images

best_fit_MI = zeros(numel(cinds), 1);
best_fit_corr = zeros(numel(cinds), 1);
for ID = cinds
    loc = max2d(MI_grid_all(:, :, ID));
    best_fit_MI(ID) = loc;
    loc = max2d(abs(corrCoeff_grid_all(:, :, ID)));
    best_fit_corr(ID) = loc;
end

%%
a = 200;
w = 3500;

im_photon = load(fullfile(folder, sprintf(fn_photons_in, 1)));

im_tm = load(fullfile(folder, sprintf(fn_tm_in, 1)));

tm = mat2gray(im_tm);
phot = mat2gray(im_photon);
corrected = spcDistortionCorrection(tm, a, w);
MI_initial = mi(tm, phot);
CC_initial = corrcoef(tm, phot);
CC_initial = CC_initial(1, 2);
MI = mi(corrected, phot);
CC = corrcoef(corrected, phot);
CC = CC(1, 2);

figure
subplot(311)
imshow(tm)
title(['Initial TM, ', 'MI : ', num2str(MI_initial), ', CC : ', num2str(CC_initial)])
subplot(312)
imshow(phot)
title('Photon')
subplot(313)
imshow(corrected)
title(['Corrected TM, ', 'MI : ', num2str(MI), ', CC : ', num2str(CC)])
% subplot(414)
% xv = 1:size(tm,2);
% plot(xv,a*sin(1/w*(xv-numel(xv)/2)))