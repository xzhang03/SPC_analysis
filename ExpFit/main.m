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
            [photCount{f}, photArrival{f}] = read_sdt(fullfile(sliceDir{dirID},currFiles.name{f}),binFactor,'matfile');
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

%% Fit individual images pooled
x=bins;

g = fittype(modelF,'independent','t');
lower = [0,6.5,1,1,0.5,1.2];
upper = [0.2,7.2,1e5,1e5,1.5,3];
opts = {'METHOD','NonlinearLeastSquares','Display','off',...
    'Robust','LAR','Lower',lower,'Upper',upper,'StartPoint',mean([lower;upper])};

for acq = 1:size(allDat.photCount,2)
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
lower = [globSrc(1),globSrc(2),0.5,0.5,globSrc(3)-0.1,globSrc(4)-0.1];
upper = [globSrc(1),globSrc(2),1e2,1e2,globSrc(3)+0.1,globSrc(4)+0.1];
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

% %% Trying to convert thingsss
fieldsArray = {'A','B','tau1','tau2'};
pixParams = obj2struct(pixFits,fieldsArray);
temp=num2cell([pixParams.A]./[pixParams.B]);
[pixParams.Ratio]=temp{:};
clear temp

fieldsArray = {'sse','rsquare','dfe','adjrsquare','rmse'};
pixRes = obj2struct(pixGOFs,fieldsArray);


figure
subplot(221)
histogram([pixParams.A])
title('A')

subplot(222)
histogram([pixParams.B])
title('B')

subplot(223)
histogram([pixParams.tau1])
title('tau1')

subplot(224)
histogram([pixParams.tau2])
title('tau2')
% %% Mean arrival time for all pixels
% for acq = 1:size(fitQ,1)
%     fprintf('# Exp %d...\n',acq)
%     Im = squeeze(fitQ(acq,:,:,:));
%     sizeIm=size(Im);
%     for idx = 1:prod(sizeIm(1:2))
%         [i,j]=ind2sub(sizeIm(1:2),idx);
%         pixParams(acq,i,j).meanT =mean(squeeze(Im(i,j,:)).*bins);
%     end
% end
% disp('done')

