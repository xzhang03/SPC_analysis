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

%% Registration
rtype = 'translation';
moving = flim;
fixed = sbx;

[optimizer,metric] = imregconfig('multimodal');
disp('## Regular registration...')
tform = imregtform(moving,fixed,rtype,optimizer, metric);
disp('    Shifts : ')
disp(tform.T(3,1:2))
disp('## Phase correlation...')
tformP = imregcorr(moving,fixed,rtype);
disp('    Shifts : ')
disp(tformP.T(3,1:2))

Rfixed = imref2d(size(fixed));
flimReg = imwarp(moving,tform,'OutputView',Rfixed);
flimRegPhase = imwarp(moving,tformP,'OutputView',Rfixed);

% %% Plot output of registration
% figure
% subplot(311)
% imshowpair(sbx,moving);
% title('Initial')
% subplot(312)
% imshowpair(sbx,flimReg);
% title('Regular registration')
% subplot(313)
% imshowpair(sbx,flimRegPhase);
% title('Phase registration')


%% GUI tests

figure
imshow(flimRegPhase)

r1 = ginput(4);
%%
flimCorrected = spcDistortionCorrection(flimRegPhase,3000,800);

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