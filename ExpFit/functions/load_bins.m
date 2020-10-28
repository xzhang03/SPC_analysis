function bins=load_bins(path)
if nargin<1
    path = "C:\Users\hfluhr\Documents\SPC_analysis\ExpFit\200305_SZ333_run1_c1_data_trace.asc";
end
bins = load(path);
bins = bins(:,1);
end

