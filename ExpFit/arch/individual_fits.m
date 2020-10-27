clearvars -except data Q

rootdir = 'H:\2p\stephen\';
filelist = dir(fullfile(rootdir,'**/*.sdt'));
filelist = struct2table(filelist);
filelist = filelist(~contains(filelist.folder,'blur'),:);
bins=load_bins();

% group slices and GRIN experiments
slices = filelist(contains(filelist.folder,'slice'),:);
% grins = filelist(contains(filelist.folder,'run'),:);
clear filelist rootdir

sliceDir = unique(slices.folder,'stable');

%% focus on 1 imaging session and get all photon counts
disp('######### Loading Data #########')

if ~exist('Q','var')
    data = struct('folder',{},'files',{},'photCount',{},'photArrival',{},'acc',{});
    for dirID = 1:3%numel(sliceDir)
        disp(['########## DIR : ',num2str(dirID),' ##########'])
        data(dirID).folder = sliceDir{dirID};
        currFiles = slices(strcmp(slices.folder,data(dirID).folder),:);
        data(dirID).files = currFiles;
        parfor f = 1:size(currFiles,1)
            [photCount{f}, photArrival{f}] = read_sdt(fullfile(sliceDir{dirID},currFiles.name{f}),10);
        end
        data(dirID).photCount = horzcat(photCount{:});
        data(dirID).photArrival = cat(4,photArrival{:});
        data(dirID).acc = sum(data(dirID).photCount,2);
        
        clear photCount photArrival
    end
    clear dirID
    Q=cat(4,data(:).photArrival);
    Q=permute(Q,[4 1 2 3]);
    data=rmfield(data,'photArrival');
end

%% Complete model definition
disp('######### Model Definition #########')

expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same');

lower = [0,5,0,0,1.2,1.5];
upper = [0.2,max(bins),1e2,1e2,1.5,2];
startP = [1,1,1,1,1,1];

g = fittype(modelF,'independent','t');

opts = {'METHOD','NonlinearLeastSquares',...
    'Display','off','Robust','off','Lower',lower,'Upper',upper,...
    'StartPoint',startP,'MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-6};

%% Fitting
disp('######### Fitting #########')


x=bins;
sizeQ = size(Q);
fitTens=cell(sizeQ(1:3));
sizePix=sizeQ(2:3);

%% Creating subsample
% Uniform pixel sampling
% N = 100;
% s = randsample(prod(sizeQ(1:3)),N);
% clear ID
% [ID.exp,ID.x,ID.y] = ind2sub(sizeQ(1:3),s);
% ID = [ID.exp,ID.x,ID.y]; 
% subQ = zeros(N,256);
% for i = 1:N
%     subQ(i,:)=squeeze(Q(ID(i,1),ID(i,2),ID(i,3),:));
% end

% Uniform images sampling
N = 10;
s = randsample(sizeQ(1),N);
subQ=Q(1,:,:,:);

%% Annoying loop - reshaped array
fitMat=cell(N);
tic
parfor idx = 1:N
    [i, j] = ind2sub(sizePix,idx);
    y = subQ(idx,:)';
    fo = fitoptions(opts{:});
    fitMat{idx} = fit(x,y,g,fo);
end
toc

%% plotting results
figure
hold on
for idx = 1:N
    y = subQ(idx,:)';
    clf
    f=fitMat{idx};
    plot(x,y,x,modelF(f.w,f.shift,f.A,f.B,f.tau1,f.tau2,x),'.')
%     p.Color(4)=1e-1;
%     legend('Signal','Fit')
    pause(1)
    drawnow
end
%% Annoying loop

for exp = 1:size(subQ,1)
    tic
%     if mod(exp,10)==0
    fprintf('# Exp %d...\n',exp)
%     end
%     Im = squeeze(subQ(exp,:,:,:));
%     fitMat=cell(sizePix);
    for idx = 1:prod(sizePix)
        [i, j] = ind2sub(sizePix,idx);
%         y = squeeze(Im(i,j,:));
%         fo = fitoptions(opts{:});
%         fitMat{idx} = fit(x,y,g,fo);
        stA(i,j)=fitMat{idx};
    end
    fitTens{exp}=fitMat;
    toc
end
disp('done')

%% plotting results
% figure
% hold on
% for idx = 1:prod(sizeQ(1:3))
%     [exp,i,j]=ind2sub(sizeQ(1:3),idx);
%     y=squeeze(Q(exp,i,j,:));
% %     f=fitTens{exp}{i,j};
% %     plot(x,y,x,modelF(f.w,f.shift,f.A,f.B,f.tau1,f.tau2,x),'.')
%     p=plot(x,y,'k');
%     p.Color(4)=1e-3;
% %     legend('Signal','Fit')
% %     drawnow
% end
