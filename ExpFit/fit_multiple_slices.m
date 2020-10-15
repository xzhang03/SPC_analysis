clearvars -except akarMice

rootdir = 'H:\2p\stephen\';
filelist = dir(fullfile(rootdir,'**/*.sdt'));
filelist = struct2table(filelist);
filelist = filelist(~contains(filelist.folder,'blur'),:);

%% group slices and GRIN experiments
slices = filelist(contains(filelist.folder,'slice'),:);
grins = filelist(contains(filelist.folder,'run'),:);

sliceDir = unique(slices.folder,'stable');
sliceDir = sliceDir(contains(sliceDir,akarMice));

% %% afas
% [phot_count,out_image] = read_sdt(fullfile(slice_dir{1},slices.name{1}));
% 
% y=squeeze(out_image(100:110,100:110,:));
% y = squeeze(sum(sum(y,1),2));
% plot(y)
%% get all photon counts
if ~exist('Q','var')
    data = struct('folder',{},'files',{},'photCount',{},'photArrival',{},'acc',{});
    for dirID = 1:numel(sliceDir)
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

%% Show photon count, decay
figure
hold on
plot([data.acc])

%% Aggregate all data
allDat.photCount = horzcat(data(:).photCount);
allDat.acc = sum(allDat.photCount,2);

plot_decay(allDat.acc);
bins=load_bins();

%% Reduced fit, ignoring IRF
x = bins(25:200);
y = allDat.acc(25:200);

singleExp = @(A,tau,s,t) A*exp(-((t-s)./(tau)));
dualExp = @(A,B,tau1,tau2,s,t) A*exp(-(t-s)./tau1)+B*exp(-(t-s)./tau2);

currentModel = dualExp;
if isequal(currentModel,singleExp)
    lower = [1e6,1,-10];
    upper = [1e9,Inf,10];
elseif isequal(currentModel,dualExp)
    lower = [1e3,1e3,1,1,-10];
    upper = [1e9,1e9,Inf,Inf,-10];
end

g = fittype(currentModel,'independent','t');
opts = fitoptions('METHOD','NonlinearLeastSquares',...
    'Display','final','Robust','Bisquare','Lower',lower,'Upper',upper,...
    'MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-6);
fitted = fit(x,y,g,opts);
disp(fitted)

figure
plot(x,y,x,fitted(x));

%% Complete model
expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,offset,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same')+offset;

lower = [0,6.5,1,1,0.5,1.2,0];
upper = [0.2,7.2,1e10,1e10,1.5,3,0];
% weights = [zeros(1,17) ones(1,209) zeros(1,30)];
weights = ones(1,256);

x=bins;
y=allDat.acc;
% y=squeeze(Q(2,25,50,:));
g = fittype(modelF,'independent','t');
opts = fitoptions('METHOD','NonlinearLeastSquares','Weights',weights,...
    'Display','final','Robust','off','Lower',lower,'Upper',upper,...
    'MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-6);
f = fit(x,y,g,opts);
disp(f)
figure
plot(x,y,x,modelF(f.w,f.shift,f.A,f.B,f.tau1,f.tau2,f.offset,x),'.')
legend('Signal','Fit')

%% NEW model
% Hfunc = @(t0,tauG,tau,t) (1./2).*exp(((tauG.^2)/(2*tau))-((t-t0)/tau)).*(1+erf((tauG.^2-tau*(t-t0))./(sqrt(2)*tau*tauG)));
% newModel = @(t0,tauG,A,B,tau1,tau2,t) (A*Hfunc(t0,tauG,tau1,t)+B*Hfunc(t0,tauG,tau2,t));
% Hfunc = @(t0,tauG,tau,t) (1./2).*exp(((tauG.^2)/(2*tau))-((t-t0)/tau)).*(1+erf((tauG.^2-tau*(t-t0))./(sqrt(2)*tau*tauG)));
% newModel = @(t0,tauG,A,B,tau1,tau2,t) (A*Hfunc(t0,tauG,tau1,t)+B*Hfunc(t0,tauG,tau2,t));

% t0 = 10;
% tauG = 0.1;
% A = 1;
% B = 1;
% tau1 = 1;
% tau2 = 2;
% figure
% subplot(121)
% plot(x,Hfunc(t0,tau2,tauG,x))
% subplot(122)
% plot(x,modelF(t0,tauG,A,B,tau1,tau2,x))
% %%
lower = [2,0,1e4,1e4,1,1.5];
upper = [2,10,1e9,1e9,2,3];
% weights = [zeros(1,17) ones(1,209) zeros(1,30)];
weights = ones(1,234);
% weights=y;

x=bins;%(23:end);
y=allDat.acc;%(23:end);
g = fittype(newModel,'independent','t');
opts = fitoptions('METHOD','NonlinearLeastSquares','Weights',weights,...
    'Display','final','Robust','LAR','Lower',lower,'Upper',upper,...
    'MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-6);
f = fit(x,y,g,opts);
disp(f)
clf
plot(x,y,x,newModel(f.t0,f.tauG,f.A,f.B,f.tau1,f.tau2,x),'.')
legend('Signal','Fit')

%% Complete model - individual acquisitions
expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,offset,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same')+offset;

lower = [0,5,1e3,1e3,1,1.5,0];
upper = [1,max(bins),1e9,1e9,1.5,2,0];
% weights = [zeros(1,17) ones(1,209) zeros(1,30)];
weights = ones(1,256);

x=bins;

for acq = 1%:numel(dat)
    disp(['################## Folder : ', num2str(acq),' ##################'])
    y = dat(acq).photCount(:,1);
    g = fittype(modelF,'independent','t');
    opts = fitoptions('METHOD','NonlinearLeastSquares','Weights',weights,...
        'Display','notify','Robust','Bisquare','Lower',lower,'Upper',upper,...
        'MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-6);
    f{acq} = fit(x,y,g,opts);
    disp(f{acq})
    clf
    plot(x,y,x,modelF(f{acq}.w,f{acq}.shift,f{acq}.A,f{acq}.B,f{acq}.tau1,f{acq}.tau2,f{acq}.offset,x),'.')
    legend('Signal','Fit')
    title(['Folder : ', num2str(acq)])
    drawnow
end

%% Single Exp - individual acquisitions
expDecay = @(A,tau,t) heaviside(t).*(A.*exp(-t./tau));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,tau,offset,t) mean(diff(t))*conv(expDecay(A,tau,t),IRF(w,shift,t),'same')+offset;

lower = [0,5,1e3,1,0];
upper = [1,max(bins),1e9,3,0];
% weights = [zeros(1,17) ones(1,209) zeros(1,30)];
weights = ones(1,256);

x=bins;

for acq = 1:numel(dat)
    disp(['################## Folder : ', num2str(acq),' ##################'])
    y = dat(acq).acc;
    g = fittype(modelF,'independent','t');
    opts = fitoptions('METHOD','NonlinearLeastSquares','Weights',weights,...
        'Display','notify','Robust','Bisquare','Lower',lower,'Upper',upper,...
        'MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-9,'TolX',1e-6);
    f{acq} = fit(x,y,g,opts);
    disp(f{acq})
    clf
    plot(x,y,x,modelF(f{acq}.w,f{acq}.shift,f{acq}.A,f{acq}.tau,f{acq}.offset,x),'.')
    legend('Signal','Fit')
    title(['Folder : ', num2str(acq)])
    drawnow
end

%% Test model fct shape
expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,offset,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same')+offset;


plotC=true;
if plotC
    w=0.05;
    shift=7.1;
    A=2.4e7;
    B=2.4e7;
    tau1=1.600;
    tau2=1.80;
    tv = [-flip(bins(1:64));bins];
    figure
    plot(tv,expDecay(A,B,tau1,tau2,tv),'.')
    hold on
    plot(tv,modelF(w,shift,A,B,tau1,tau2,0,tv))
    plot(tv,normpdf(tv,0,w))
    plot(bins,allDat.acc)
    legend('Ideal Exp decay','Modeled signal','modeled IRF','Pooled data')
end