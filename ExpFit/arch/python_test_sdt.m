clearvars -except akarMice

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

%%
path = fullfile(slices.folder{1},slices.name{20});

%%
if count(py.sys.path,'')==0
    insert(py.sys.path,int32(0),'');
end