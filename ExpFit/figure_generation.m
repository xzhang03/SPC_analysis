% Script to get figures from the main.m script in same directory, most
% figures require to have some variables computed in that other script


%% Global fit figure
% figure
% plot(x,y,'LineWidth',1)
% hold on
% plot(x,modelF(f.w,f.shift,f.A,f.B,f.tau1,f.tau2,x),'r.')
% legend('Signal','Fit')
% xlabel('[ns]')
% ylabel('Photon Count')
% title('All pixels pooled together')

%% figure of single pixel fit
% i=40;
% j=40;
% 
% currFit=pixParams(1,i,j);
% figure
% plot(bins,squeeze(fitQ(1,i,j,:)))
% hold on
% plot(bins,modelF(globSrc(1),globSrc(2),currFit.A,currFit.B,globSrc(3),globSrc(4),bins),'LineWidth',2)
% xlabel('[ns]')
% ylabel('Photon Count')
% title('Single Pixel Fit')


%% Look at Ratio over images
% currentImg = 1;
% 
% figure
% subplot(2,1,1)
% plotStructField(pixParams(currentImg,:,:),'Ratio')
% title('Ratio')
% subplot(2,1,2)
% plotStructField(pixRes(currentImg,:,:),'rmse')
% title('RMSE')

% figure
% heatmap(getFieldArray(pixParams(currentImg,:,:),'Ratio'));

%% Test model fct shape
expDecay = @(A,B,tau1,tau2,t) heaviside(t).*((A.*exp(-t./tau1)+B.*exp(-t./tau2)));
IRF = @(w,shift,t) normpdf(t-shift,0,w);
modelF = @(w,shift,A,B,tau1,tau2,t) mean(diff(t))*conv(expDecay(A,B,tau1,tau2,t),IRF(w,shift,t),'same');

w=0.1;
shift=7.1;
A=1;
B=1;
tau1=0.7;
tau2=2;
tv = [-flip(bins(1:64));bins];
figure
subplot(311)
plot(tv,expDecay(A,B,tau1,tau2,tv),'g.')%,'LineWidth',2)
legend('Exponential Decay')
subplot(312)  
plot(tv,normpdf(tv,0,w),'LineWidth',2)
legend('Instrument Response Function')
subplot(313)
plot(tv,modelF(w,shift,A,B,tau1,tau2,tv),'r','LineWidth',2)
xlabel('[ns]')
legend('Modeled signal')
%     plot(bins,allDat.acc)

%% Figure sum of two exponentials
% 
% figure
% plot(bins,exp(-bins./tau1),'--','LineWidth',1)
% hold on
% plot(bins,exp(-bins./tau2),'--','LineWidth',1)
% plot(bins,0.5*exp(-bins./tau1)+0.5*exp(-bins./tau2),'LineWidth',2)
% xlabel('[ns]')
% ylabel('Photon Count')
% legend({'No Interaction','Interaction','Total Signal'})
% xlim([-1 14])

%% scatter plot
% figure
% scatter([pixParams.meanT],[pixParams.Ratio])
% xlabel('Mean Arrival Time')
% ylabel('Ratio')
% title('Ratio vs Mean Arrival Time for few frames')