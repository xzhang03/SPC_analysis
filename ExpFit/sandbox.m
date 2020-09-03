clear
tm_path = "H:\2p\stephen\SZ334\200624_SZ334\FLIM\200624_SZ334_run1\200624_SZ334_run1_c1_t1.asc";
sdt_path = 'H:\2p\stephen\SZ334\200624_SZ334\FLIM\200624_SZ334_run1\200624_SZ334_run1_c1.sdt';

im = load(tm_path);
sdt = bfopen(sdt_path);

%% Manip sdt file
out = cat(3,sdt{1}{:,1});

phot_count = squeeze(sum(sum(out,1),2));

%%
figure
plot(phot_count)
xlabel('time bin')
ylabel('nb of photons')