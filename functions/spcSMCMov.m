function spcSMCMov()
fp = '\\nasquatch\data\2p\stephen\SZ1086\230603_SZ1086\FLIM\230603_SZ1086_run1\';
fns = dir(fullfile(fp, '*smc.mat'));

for i = 1 : 90
    smc = load(fullfile(fns(i).folder, fullfile(fns(i).name)));
    smc = smc.smc;
    smc = spcSMCcrop(smc, [360, 20, 950, 500]);
    f = spcSMCphotonframe(smc);
    f = binxy(f, 2);
    
    if i == 1
        stack = double(repmat(f, [1 1 90]));
    else
        stack(:,:,i) = f;
    end
    
end

writetiff(stack, fullfile(fp, 'smc_demo', '230603_SZ1086_run1_photons.tif'), 'double');

end