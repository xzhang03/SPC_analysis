clearvars -except akarMice data Q

rootdir = 'H:\2p\stephen\';
filelist = dir(fullfile(rootdir,'**/*.sdt'));
filelist = struct2table(filelist);
filelist = filelist(~contains(filelist.folder,'blur'),:);

%% group slices and GRIN experiments
slices = filelist(contains(filelist.folder,'slice'),:);
grins = filelist(contains(filelist.folder,'run'),:);

sliceDir = unique(slices.folder,'stable');
if~exist('akarMice','var') load('akarMice.mat'); end
sliceDir = sliceDir(contains(sliceDir,akarMice));

%% get all photon counts
binFactor = 5;
frames2drop = 20;
sliceDir=sliceDir([1:3,5:6,8:11]);
if ~exist('Q','var')
    data = struct('folder',{},'files',{},'photCount',{},'photArrival',{},'acc',{});
    for dirID = 1:numel(sliceDir)
        disp(['########## DIR : ',num2str(dirID),' ##########'])
        data(dirID).folder = sliceDir{dirID};
        currFiles = slices(strcmp(slices.folder,data(dirID).folder),:);
        currFiles = currFiles(frames2drop:end,:);
        data(dirID).files = currFiles;
        parfor f = 1:size(currFiles,1)
            [photCount{f}, photArrival{f}] = read_sdt(fullfile(sliceDir{dirID},currFiles.name{f}),binFactor);
        end
        data(dirID).photCount = horzcat(photCount{:});
        data(dirID).photArrival = cat(4,photArrival{:});
        data(dirID).acc = sum(data(dirID).photCount,2);
        
        clear photCount photArrival
    end
    clear dirID
    Q=cat(4,data(:).photArrival);
    data=rmfield(data,'photArrival');
    Q=permute(Q,[4 1 2 3]);
end

%% Aggregate all data
allDat.photCount = horzcat(data(:).photCount);
allDat.acc = sum(allDat.photCount,2);


%%
bins=load_bins();
x=bins;


%% Complete model - Fit all aggregated data
expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same');

lower = [0,6.5,1,1,0.5,1.2];
upper = [0.2,7.2,1e10,1e10,1.5,3];
opts = {'METHOD','NonlinearLeastSquares','Display','off',...
    'Robust','LAR','Lower',lower,'Upper',upper,'StartPoint',mean([lower;upper])};

y=allDat.acc;
g = fittype(modelF,'independent','t');
fo = fitoptions(opts{:});
[f,gof] = fit(x,y,g,fo);
disp(f)
%%
figure
plot(x,y,'LineWidth',1)
hold on
plot(x,modelF(f.w,f.shift,f.A,f.B,f.tau1,f.tau2,x),'r.')
legend('Signal','Fit')
xlabel('[ns]')
ylabel('Photon Count')
title('All pixels pooled together')

%% Fit individual images pooled
x=bins;

g = fittype(modelF,'independent','t');
lower = [0,6.5,1,1,0.5,1.2];
upper = [0.2,7.2,1e5,1e5,1.5,3];
opts = {'METHOD','NonlinearLeastSquares','Display','off',...
    'Robust','LAR','Lower',lower,'Upper',upper,'StartPoint',mean([lower;upper])};

parfor acq = 1:size(allDat.photCount,2)
    y = allDat.photCount(:,acq);
    fo = fitoptions(opts{:});
    [AcqFits{acq}, AcqGOF{acq}] = fit(x,y,g,fo);
end

%% Parsing data from previous loop and comparing results for global parameters
fieldsArray = {'w','shift','A','B','tau1','tau2'};
imgFits = obj2struct(AcqFits,fieldsArray);
fieldsArray = {'sse','rsquare','dfe','adjrsquare','rmse'};
imgGOFs = obj2struct(AcqGOF,fieldsArray);

Parameters = {'tauG';'t0';'tau1';'tau2'};
SingleFit = [f.w;f.shift;f.tau1;f.tau2];
IndivFits = mean([imgFits.w;imgFits.shift;imgFits.tau1;imgFits.tau2],2);
globParamsT = table(Parameters,SingleFit,IndivFits);
disp(globParamsT)
clear fieldsArray Parameters

%% Fit individual pixels
disp('######### Fitting all Pixels #########')
x=bins;
fitQ = Q(1,:,:,:);
sizeQ = size(fitQ);
pixFits=cell(sizeQ(1:3)); pixGOFs=cell(sizeQ(1:3));
sizePix=sizeQ(2:3);

globSrc = SingleFit;
lower = [globSrc(1),globSrc(2),1,1,globSrc(3),globSrc(4)];
upper = [globSrc(1),globSrc(2),1e5,1e5,globSrc(3),globSrc(4)];
opts = {'METHOD','NonlinearLeastSquares','Display','off',...
    'Robust','LAR','Lower',lower,'Upper',upper,'StartPoint',mean([lower;upper])};

for acq = 1:size(fitQ,1)
    fprintf('# Exp %d...\n',acq)
    Im = squeeze(fitQ(acq,:,:,:));
    fitMat=cell(sizePix); fitGOFs=cell(sizePix);
    parfor idx = 1:prod(sizePix)
        [i, j] = ind2sub(sizePix,idx);
        y = squeeze(Im(i,j,:));
        fo = fitoptions(opts{:});
        [fitMat{idx},fitGOFs{idx}] = fit(x,y,g,fo);
    end
    pixFits(acq,:,:)=fitMat; pixGOFs(acq,:,:)=fitGOFs;
end
disp('done')

%% Trying to convert thingsss
fieldsArray = {'A','B'};
pixParams = obj2struct(pixFits,fieldsArray);
temp=num2cell([pixParams.A]./[pixParams.B]);
[pixParams.Ratio]=temp{:};
clear temp

fieldsArray = {'sse','rsquare','dfe','adjrsquare','rmse'};
pixRes = obj2struct(pixGOFs,fieldsArray);
return

%% figure of single pixel fit
i=40;
j=40;

currFit=pixParams(1,i,j);
figure
plot(bins,squeeze(fitQ(1,i,j,:)))
hold on
plot(bins,modelF(globSrc(1),globSrc(2),currFit.A,currFit.B,globSrc(3),globSrc(4),bins),'LineWidth',2)
xlabel('[ns]')
ylabel('Photon Count')
title('Single Pixel Fit')


%% Look at Ratio over images
currentImg = 1;

figure
% subplot(2,1,1)
plotStructField(pixParams(currentImg,:,:),'Ratio')
title('Ratio')
% subplot(2,1,2)
% plotStructField(pixRes(currentImg,:,:),'rmse')
% title('RMSE')

%%
figure
heatmap(getFieldArray(pixParams(currentImg,:,:),'Ratio'));

%% Mean arrival time for all pixels
for acq = 1:size(fitQ,1)
    fprintf('# Exp %d...\n',acq)
    Im = squeeze(fitQ(acq,:,:,:));
    sizeIm=size(Im);
    for idx = 1:prod(sizeIm(1:2))
        [i,j]=ind2sub(sizeIm(1:2),idx);
        pixParams(acq,i,j).meanT =mean(squeeze(Im(i,j,:)).*bins);
    end
end
disp('done')

%% scatter plot
figure
scatter([pixParams.meanT],[pixParams.Ratio])
xlabel('Mean Arrival Time')
ylabel('Ratio')
title('Ratio vs Mean Arrival Time for few frames')

%% view
dirID=1; f=50;
currFiles = slices(strcmp(slices.folder,data(dirID).folder),:);
path =fullfile(sliceDir{dirID},currFiles.name{f});
 
% grinDir = unique(grins.folder,'stable');
% grinDir = grinDir(contains(grinDir,akarMice));
% currFiles=grins(strcmp(grins.folder,grinDir(1)),:);
% path = fullfile(currFiles.folder{1},currFiles.name{1});
%%
% [photCount, photArrival] = read_sdt(path,binFactor);
% sdt = bfopen(path);
% outImage = cat(3,sdt{1}{:,1});
% 
% % try to get mean arrival time
% A=double(squeeze(outImage));
% sizeA=size(A);
% meanT = zeros(sizeA(1:2));
% parfor idx = 1:prod(sizeA(1:2))
%     [i,j]=ind2sub(sizeA(1:2),idx);
%     meanT(idx) = mean(squeeze(A(i,j,:)).*bins);
% end
% %
% figure
% imshow(mat2gray(meanT))

%% Trying to reshape the sdt output file
out=reshape([meanT(2:2:end,:) meanT(1:2:end,end:-1:1)],size(meanT,1),size(meanT,2));
figure
imshow(mat2gray(meanT))


% %% Test model fct shape
% expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
% IRF = @(w,shift,t) normpdf(t-shift,0,w);
% modelF = @(w,shift,A,B,tau1,tau2,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same');
% 
% 
% plotC=true;
% if plotC
%     w=0.1;
%     shift=7.1;
%     A=1;
%     B=1;
%     tau1=0.7;
%     tau2=2;
%     tv = [-flip(bins(1:64));bins];
%     figure
%     subplot(311)
%     plot(tv,expDecay(A,B,tau1,tau2,tv),'g.')%,'LineWidth',2)
%     legend('Exponential Decay')
%     subplot(312)  
%     plot(tv,normpdf(tv,0,w),'LineWidth',2)
%     legend('Instrument Response Function')
%     subplot(313)
%     plot(tv,modelF(w,shift,A,B,tau1,tau2,tv),'r','LineWidth',2)
%     xlabel('[ns]')
%     legend('Modeled signal')
% %     plot(bins,allDat.acc)
% end
% %%
% 
% figure
% plot(bins,exp(-bins./tau1),'--','LineWidth',1)
% hold on
% plot(bins,exp(-bins./tau2),'--','LineWidth',1)
% plot(bins,0.5*exp(-bins./tau1)+0.5*exp(-bins./tau2),'LineWidth',2)
% xlabel('[ns]')
% ylabel('Photon Count')
% legend({'No Interaction','Interaction','Total Signal'})
xlim([-1 14])