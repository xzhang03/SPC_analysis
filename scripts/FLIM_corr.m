%% Look
data = stkread2('F:\2p\stephen\SZ333\200609_SZ333\FLIM\Output2\C1-AVG_200609_SZ333_run1_odd_even cropped.tif');
data2 = reshape(data, [242*347, 2]);
ft2(data2(:,1),10);

%% Low pass
dl = designfilt('lowpassfir', 'StopbandFrequency', Fpass, 'PassbandFrequency', ...
    Fstop, 'StopbandAttenuation', Astop, 'PassbandRipple', Apass, 'SampleRate',...
    Fs, 'DesignMethod', 'equiripple');

fvtool(dl)

data3_1d = filter(dl, data2(:,1));
data3_2d = filter(dl, data2(:,2));
imshow(reshape(data3_1d, [242, 347]), [])

%% High pass
Fstop = 0.08;
Fpass = 0.1;
Apass = 0.5;
Astop = 60;
Fs = 10;

d = designfilt('highpassfir', 'StopbandFrequency', Fstop, 'PassbandFrequency', ...
    Fpass, 'StopbandAttenuation', Astop, 'PassbandRipple', Apass, 'SampleRate',...
    Fs, 'DesignMethod', 'equiripple');

fvtool(d)
data3_1 = filter(d, data2(:,1));
data3_2 = filter(d, data2(:,2));
imshow(reshape(data3_1, [242, 347]), [])

%% High pass samples
hp1 = mat2gray([reshape(data3_1, [242, 347]), reshape(data3_2, [242, 347])]);
hp2 = mat2gray([reshape(data3_1d, [242, 347]), reshape(data3_2d, [242, 347])]);
imtool([hp1;hp2-0.9],[])


%% Correlation
chp = corr(data3_1(:),data3_2(:));

%% Gaussian pass
data4_1 = data(:,:,1) - imgaussfilt(data(:,:,1), 12);
data4_2 = data(:,:,2) - imgaussfilt(data(:,:,2), 12);
imtool(data4_1, [])

%% Gaussian samples
gs = [mat2gray([data4_1,data4_2]);
    mat2gray([imgaussfilt(data(:,:,1),12), imgaussfilt(data(:,:,2),12)])];
imtool(gs,[])

%% Correlation
cgf = corr(data4_1(:),data4_2(:));

%% Shuffle
nitr = 10000;
maxind = 242*347;
vs = zeros(nitr,2);

hbar = waitbar(0, 'shuffling');
for i = 1 : nitr
    if mod(i,1000) == 0
        waitbar(i/nitr);
    end
    
    inds = randi(maxind, [maxind,1]);
    vs(i,1) = corr(data3_1(:),data3_2(inds));
    vs(i,2) = corr(data4_1(:),data4_2(inds));
end
close(hbar)

p_hp = mean(vs(:,1) > chp);
p_gf = mean(vs(:,2) > cgf);