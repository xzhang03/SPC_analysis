
pathSDT = 'H:\2p\stephen\SZ309\FLIM\200128_SZ309\200128_SZ309_slice1\200128_SZ309_slice1_c050.sdt';
pathIMG = 'H:\2p\stephen\SZ309\FLIM\200128_SZ309\200128_SZ309_slice1\200128_SZ309_slice1_c050.img';
pathTiff = "H:\2p\stephen\SZ309\FLIM\200128_SZ309\200128_SZ309_slice1\200128_SZ309_slice1_c020_intensity_image.tif";
%%
sdtO = bfopen(pathSDT);
outImage = cat(3,sdtO{1}{:,1});
sdt = squeeze(sum(outImage,3));

%%
tf = Tiff(pathTiff);
img = double(read(tf));

%% Histogram of pixel values to find correspondance
figure
subplot(211)
histogram(sdt,10)
lim=get(gca,'YLim');
subplot(212)
histogram(img,10)
set(gca,'YLim',lim);

%%
disp(max2d(sdt))
disp(max2d(img))

%%
[C,ia,ic]=unique(img);
a_counts=accumarray(ic,1);
value_countsIMG=[C,a_counts];
[C,ia,ic]=unique(sdt);
a_counts=accumarray(ic,1);
value_countsSDT=[C,a_counts];
value_countsSDT=[value_countsSDT;44,0];

%%
vT=table(value_countsSDT(:,1),value_countsSDT(:,2), value_countsIMG(:,2),'VariableNames',{'Value';'SDT';'IMG'})