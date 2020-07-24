clear
mouse = 'SZ334';
date = 200624;
run = 1;
user = 'stephen';

%% Get paths
spcpaths = spcPath(mouse, date, run, 'user',user);
inds = spcpaths.cinds;

%% Loading photon and tm files
im_photon = load(fullfile(spcpaths.fp, sprintf(spcpaths.photons_in,ind)));

im_tm = load(fullfile(spcpaths.fp, sprintf(spcpaths.tm_in,ind)));

%%
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

%% cosine fitting
w_array = ceil(size(tm,2))/2:200:200000;
a_array = 1:200:100000;
[A,W] = ndgrid(a_array,w_array);


MI_grid = zeros(length(w_array),length(a_array));
corrCoeff_grid = zeros(length(w_array),length(a_array));
binxy
parfor IDX = 1:numel(A)
    a = A(IDX);
    w = W(IDX);
    corrected = sine_unwarp(tm,a,w);
    MI_grid(IDX) = mi(corrected, phot);
    img_corr = corrcoef(corrected,phot);
    corrCoeff_grid(IDX) = img_corr(1,2);
end

%%
outgrid = mat2gray(MI_grid);
outgrid = imgaussfilt(outgrid,20);
figure
imshow(outgrid')
