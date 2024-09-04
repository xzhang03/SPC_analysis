# spc_analysis
 FLIM data analysis in MATLAB, written for Becker & Hickl FLIM system.
 
 Written by Stephen Zhang and Hugo Fluhr

### Purpose
The main purpose of this repository is to provide a package to load and process FLIM data

### Main functions
1. Load sdt files and compress them to a new format called sparse-matrix compression
2. Register FLIM files using rigid and non-rigid algorithms
3. Prepare FLIM images for segmentation
4. Apply segments to extract intensity and lifetime traces
5. Provide multiple algorithms of lifetime calculations: fitting, Tmean, IEM
6. Save tiff files of photon counts, Tmeans, IEMs
7. Correct for image-warping that is resulted from resonant scanning
8. Base analysis for FLIP data

### Folders

1. DistortionCorrection: Code to correct for image warping that is resulted from resonant scanning
2. ExpFit: Code for GLM-based exponential fitting algorithm (not commonly used)
3. Sample_data: Sample data to demonstrate the SMC (sparse-matrix compression algorithm)
4. bfmatlab: bioformat support for matlab - for loading FLIP data
5. functions: Main functions for FLIM analysis
