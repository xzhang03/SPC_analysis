## 1. SMC function list
SMC is a FLIM data algorithm developed by Stephen to save a 3D/4D data matrix in a sparse matrix form. Compression rate is typically -98-99% depending on data density. Please see sample data folder in the main directory to see example compression and decompression results

#### Apply SMC compression to 3d/4d FLIM data matrix:
```Matlab
smstruct = spcCompress(mov, varargin)
```

#### Decompress SMC files to 3d/4d FLIM data matrix:
```Matlab
mov = spcDecompress(smstruct, varargin)
```

#### Apply cellpose-generated ROI to SMC-compressed FLIM data and extract traces:
```Matlab
spcApplyCellposeROI_smc(mouse, date, run, varargin)
```

#### Demonsreg for SMC files:
```Matlab
spcSMCDemreg(mouse, date, run, varargin)
```

#### Load SMC files:
```Matlab
[output, n, smccrop] = spcSMCLoad(mouse, date, run, varargin)
```

#### XY rigid registration for SMC files:
```Matlab
spcSMCReg(mouse, date, run, varargin)
```

#### Quick analysis of SMC data:
```Matlab
output_struct = spcSMCanalysis(smc, varargin)
```

#### Bin SMC files:
```Matlab
smc = spcSMCbin(smc, xybin, tbin)
```

#### Crop SMC files:
```Matlab
smc = spcSMCcrop(smc, crop)
```

#### Generate a photon-intensity frame from SMC data:
```Matlab
frames = spcSMCphotonframe(smc, varargin)
```

## 2. Function list:

#### Make neuropil ring based on defined inner and outer diameters (in pixels):
```Matlab
[np, np_size] = MakeNPRing(im, inner, outter, im_neg, im_allow, minpixels)
```

#### A simple function to get user-drawn polygons, used to define the contour of GRIN lens:
```Matlab
[ outputim ] = getpoly( inputim , text)
```

#### Use 7zip to compress data files for long-term storage:
```Matlab
spc7zip(mice, ziplevel, varargin)
```

#### Performs area-under-curve analysis for data:
```Matlab
Areas = spcAUC(data)
```

#### Apply cellpose-generated ROI to FLIM data and extract traces:
```Matlab
spcApplyCellposeROI(mouse, date, run, varargin)
```

#### Applies ROIs to FLIM images. The slice version generates a matrix that is compatible with other slice processing (obsolete):
```Matlab
spcApplyROIslice(mouse, date, run, varargin)
```

#### Apply Sbx ROIs to FLIM images:
```Matlab
spcApplySbxROI(mouse, date, run, varargin)
```

#### Prepare an image for cellpose-based segmentation:
```Matlab
spcCellposePrep(mouse, date, run, varargin)
```

#### Automatically crops an image:
```Matlab
spcAutocrop(im, varargin)
```

#### Check xy bins for image and for trace extraction in preparation for lifetime calculations:
```Matlab
spcCheckBinxy(mouse, date, run, varargin)
```

#### Check expression level using intensity data:
```Matlab
spcCheckExpression(mouse, date, run, varargin)
```

#### Copy data to dropbox:
```Matlab
spcDBCopy(mouse, date, varargin)
```

#### Load raw sdt files (not using bioformat toolbox):
```Matlab
[mov, params] = spcLoadsdt(fpath, varargin)
```

#### Load raw sdt files (using bioformat toolbox). This code is now used for FLIP data only:
```Matlab
[mov, params] = spcLoadsdt_bf(fpath, rowinc, movdim)
```

#### Get decay trace from a 3d FLIM data file:
```Matlab
t = spcMovtrace(mov, r, c, binsz)
```

#### Calculate neuropil lifetimes:
```Matlab
spcNeuropil(mouse, date, varargin)
```

#### Paths for FLIM files:
```Matlab
spcpaths = spcPath(mouse, date, run, varargin)
```

#### Use a UI to get a quick peek of FLIM data series (sdt files). This is usually the first code to run:
```Matlab
savestruct = spcPeeksdt(mouse, date, run, varargin)
```

#### Manual ROI match across days/sessions:
```Matlab
spcROIMatchingManual(mouse, date, runs, varargin)
```

#### Manual ROI matching using a dedicated reference set:
```Matlab
spcROIMatchingManualRef(mouse, date, runs, varargin)
```

#### Manual ROI matching between sbx and FLIM data:
```Matlab
spcROIMatchingManualSbx(mousecell_spc, datecell_spc, runs_spc, mousecell_sbx, datecell_sbx, runs_sbx, optotunes_sbx, varargin)
```

#### Correct for warping that is resulted from resonant scanning:
```Matlab
spcTiffWarp(mouse, date, varargin)
```

#### Convert sdt files to tiff files directly:
```Matlab
spcTiff_sdt(mouse, date, run, varargin)
```

## 3. Obsolete functions

#### Apply ROI to pre-calculated photons/tm tiff stack (obsolete):
```Matlab
spcApplyROI(mouse, date, varargin)
```

#### FOV-wide FLIM analysis (obsolete):
```Matlab
fovstruct = spcFOVoverview(mousecell, datecell, runcell, varargin)
```

#### Get whole-FOV FLIM traces from pre-calculated tiff tiles (obsolete):
```Matlab
[tm, trace, photontrace] = spcGrossTm(mouse, date, run, varargin)
```

#### Morphological segmentation of FLIM images (obsolete):
```Matlab
spcMorphSegmentation(mouse, date, run, varargin)
```

#### Automated ROI matching across days/sessions (obsolete):
```Matlab
spcROIMatching(mouse, date, runs, varargin)
```

#### Convert SPC output to image files (obsolete):
```Matlab
spcTiff(mouse, date, run, varargin)
```

#### Uses demonsreg on the spc photon data and apply the shifts to the tm data (obsolete):
```Matlab
spcTiffDemonsReg(mouse, date, run, varargin)
```

#### XY register tiff of photon images, and apply the shifts to tm (obsolete):
```Matlab
spcTiffReg(mouse, date, run, varargin)
```
