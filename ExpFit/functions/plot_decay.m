function plot_decay(photCount)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

bins = load_bins();

ticks = 1:32:256;
figure
plot(photCount)
xlabel('time [ns]')
xticks(ticks)
xticklabels(round(bins(ticks),2))
ylabel('Photon count')
title('Fluorescence lifetime')
end

