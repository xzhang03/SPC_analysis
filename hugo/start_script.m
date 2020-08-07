machine = getenv('COMPUTERNAME');

if ~strcmp(machine,'CAMILLO')
    addpath(genpath('C:/Users/Hugo/Documents/MATLAB/SPC_analysis'))
    rmpath('C:/Users/Hugo/Documents/MATLAB/SPC_analysis/hugo/archive/')
    
    sbx_path = "C:\Users\Hugo\Documents\MATLAB\SPC_analysis\hugo\sample_images\sbx_sample.tif";
    flim_path = "C:\Users\Hugo\Documents\MATLAB\SPC_analysis\hugo\sample_images\flim_sample.tif";
else
    addpath(genpath('C:\Users\hfluhr\Documents\SPC_analysis'))
    rmpath('C:\Users\hfluhr\Documents\SPC_analysis/hugo/archive/')
    
    sbx_path = "C:\Users\hfluhr\Documents\SPC_analysis\hugo\sample_images\sbx_sample.tif";
    flim_path = "C:\Users\hfluhr\Documents\SPC_analysis\hugo\sample_images\flim_sample.tif";
end

