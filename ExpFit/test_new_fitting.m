%% Test model fct shape
expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same');
singleExp = @(tau,t) heaviside(t).*(exp(-t./tau));
% modelAlt = @(w,shift,A,B,tau1,tau2,t) mean(diff(t))*(conv(singleExp(A,tau1,t),IRF(w,shift,t),'same') +...
%     conv(singleExp(B,tau2,t),IRF(w,shift,t),'same'));

comp = @(w,shift,tau,t) mean(diff(t))*(conv(singleExp(tau,t),IRF(w,shift,t),'same'));
bins=load_bins();
 
%%
% max(modelAlt(w,shift,A,B,tau1,tau2,tv)-A*modelAlt(w,shift,1,B,tau1,tau2,tv))

w=0.0627;
shift=7.2;
tau1=0.74;
tau2=2.04;

comp1=comp(w,shift,tau1,bins);
comp2=comp(w,shift,tau2,bins);


%% test data
if ~exist('fitQ','var') 
    load fitTestData.mat
end
pixels=reshape(fitQ,[],256);
highPix =pixels(photCount>100,:);

%%
clear b dev stats
tic
for i = 1:size(highPix,1)% 1:50%size(pixels,1)
    testPix=highPix(i,:);
    [b{i},dev{i},stats{i}] = glmfit([comp1 comp2],testPix,[],'constant','off');
%     [b,dev,stats] = glmfit([comp1 comp2],testPix,[],'constant','off');
end
toc
% 
% figure
% plot(bins,testPix,bins,b(1)*comp1+b(2)*comp2)
b = cell2mat(b);

%% Stuff
% plot(bins,testPix,bins,b(1)*comp1+b(2)*comp2)
figure
subplot(131)
histogram(b(1,:),10)
title('A')
subplot(132)
histogram(b(2,:),10)
title('B')
subplot(133)
histogram(b(1,:)./b(2,:),10)
title('Ratio')
%% compare with old fitting, super slow

lower = [w,shift,0.1,0.1,tau1,tau2];
upper = [w,shift,50,50,tau1,tau2];
opts = {'METHOD','NonlinearLeastSquares','Display','final',...
    'Robust','LAR','Lower',lower,'Upper',upper,'TolX',1e-2};
x=bins;
g = fittype(modelF,'independent','t');
fo = fitoptions(opts{:});
clear f gof
tic
for i = 1:size(highPix,1)%1:size(pixels,1)
    testPix=highPix(i,:);
%     [f,gof] = fit(x,testPix',g,fo);
    [f{i},gof{i}] = fit(x,testPix',g,fo);
end
toc
% figure
% plot(bins,testPix,bins,modelF(w,shift,f.A,f.B,tau1,tau2,bins))

%% Trying to convert thingsss
fieldsArray = {'A','B','tau1','tau2'};
pixParams = obj2struct(f,fieldsArray);
temp=num2cell([pixParams.A]./[pixParams.B]);
[pixParams.Ratio]=temp{:};
clear temp

fieldsArray = {'sse','rsquare','dfe','adjrsquare','rmse'};
pixRes = obj2struct(gof,fieldsArray);


figure
subplot(131)
histogram([pixParams.A],10)
title('A')
subplot(132)
histogram([pixParams.B],10)
title('B')
subplot(133)
histogram([pixParams.Ratio],10)
title('Ratio')