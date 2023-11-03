function spcCheckBin(im)
%spcCheckBinsmc Checks bins

im2show = sum(im,3);
figure;
imshow(im2show, []);
[x, y] = ginput(1);
close(gcf);


end

