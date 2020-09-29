function bins = plot_decay(photCount)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

binRef = "C:\Users\hfluhr\Documents\SPC_analysis\ExpFit\200305_SZ333_run1_c1_data_trace.asc";
bins = load(binRef);
bins = bins(:,1);

ticks = 1:32:256;
figure
plot(photCount)
xlabel('time [ns]')
xticks(ticks)
xticklabels(round(bins(ticks),2))
ylabel('Photon count')
title('Fluorescence lifetime')
end

