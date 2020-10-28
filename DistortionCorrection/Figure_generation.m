
% Distortion Correction (Slice or GRIN depending on which script is run
%% prior to correction
figure
subplot(211)
imshow(im_sbx)
title('SBX Intensity')
subplot(212)
imshow(im_flim)
title('FLIM Intensity')

%% After correction
figure
subplot(211)
imshow(im_flim)
title('FLIM Uncorrected')
subplot(212)
imshow(flimCorrected)
title('FLIM Corrected')

%% Deformation field (while in break point xRegCorrect
figure
plot(field,'r','LineWidth',1.5)
title('Deformation Field')
xlabel('x pixels')
ylabel('y deformation')