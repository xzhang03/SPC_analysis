clear
sdt_path = 'H:\2p\stephen\SZ334\200624_SZ334\FLIM\200624_SZ334_run1\200624_SZ334_run1_c1.sdt';

[phot_count, out] = read_sdt(sdt_path);
%%
figure
plot(phot_count)
xlabel('time bin')
ylabel('nb of photons')

%% all sdt files
clear
rootdir = 'H:\2p\stephen\';
filelist = dir(fullfile(rootdir,'**/*.sdt'));
filelist = struct2table(filelist);
filelist = filelist(~contains(filelist.folder,'blur'),:);

%% group slices and GRIN experiments
slices = filelist(contains(filelist.folder,'slice'),:);
grins = filelist(contains(filelist.folder,'run'),:);

slice_dir = unique(slices.folder,'stable');


%% focus on 1 imaging session and get all photon counts
currDir = slice_dir{1};
currFiles = slices(strcmp(slices.folder,currDir),:);

parfor f = 1:size(currFiles,1)
    photCount{f} = read_sdt(fullfile(currDir,currFiles.name{f}));
end

%% 
acc = sum(cat(2,photCount{:}),2);
figure
plot(acc)
xlabel('time bin')
ylabel('nb of photons')