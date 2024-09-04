## SMC
SMC is a FLIM data algorithm developed by Stephen to save a 3D/4D data matrix in a sparse matrix form. Compression rate is typically -98-99% depending on data density. Please see sample data folder in the main directory to see example compression and decompression results

## Function list:

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
function Areas = spcAUC(data)
```

#### Apply cellpose-generated ROI to FLIM data and extract traces:
```Matlab
spcApplyCellposeROI(mouse, date, run, varargin)
```

#### Apply cellpose-generated ROI to SMC-compressed FLIM data and extract traces:
```Matlab
spcApplyCellposeROI_smc(mouse, date, run, varargin)
```

#### Apply ROI to pre-calculated photons/tm tiff stack (obsolete):
```Matlab
spcApplyROI(mouse, date, varargin)
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
function spcCheckBinxy(mouse, date, run, varargin)
```

#### Check expression level using intensity data:
```Matlab
spcCheckExpression(mouse, date, run, varargin)
```

#### Apply SMC compression to 3d/4d FLIM data matrix:
```Matlab
smstruct = spcCompress(mov, varargin)
```

#### Copy data to dropbox:
```Matlab
spcDBCopy(mouse, date, varargin)
```

#### Decompress SMC files to 3d/4d FLIM data matrix
```Matlab
mov = spcDecompress(smstruct, varargin)
```

#### FOV-wide FLIM analysis (obsolete)
```Matlab
fovstruct = spcFOVoverview(mousecell, datecell, runcell, varargin)
```

#### Get whole-FOV FLIM traces from pre-calculated tiff tiles (obsolete)
```Matlab
[tm, trace, photontrace] = spcGrossTm(mouse, date, run, varargin)
```
