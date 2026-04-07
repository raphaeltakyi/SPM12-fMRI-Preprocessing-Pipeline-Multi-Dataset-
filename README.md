# SPM12-fMRI-Preprocessing-Pipeline-Multi-Dataset-
This repository contains a MATLAB script for automated preprocessing of fMRI data using SPM12. It is designed to handle multiple datasets with different acquisition parameters, while remaining robust and flexible.

Overview
This pipeline performs standard fMRI preprocessing steps:
1. Realignment & Unwarping
2. Slice Timing Correction
3. Coregistration (T1 → Functional)
4. Segmentation (GM, WM, CSF + deformation fields)
5. Normalisation to MNI space
6. Spatial Smoothing
   
The script supports:
1. Multiple datasets with different TRs and naming conventions
2. Automatic slice number detection from NIfTI headers
3. Session-wise preprocessing (handles variability across runs)
4. Logging of preprocessing steps and errors

Folder Structure (Expected)
dataset/
├── sub-001/
│   ├── ses-01/
│   │   ├── anat/
│   │   │   └── *T1w.nii
│   │   └── func/
│   │       └── *bold*.nii
│   ├── ses-02/
│   │   └── ...


Configuration
Update the following paths in the script:
dataset1_dir = 'path/to/dataset1';
dataset2_dir = 'path/to/dataset2';
output_dir   = 'path/to/output';

TR_dataset1 = XX;   % e.g., 2.0
TR_dataset2 = XX;

Other adjustable parameters:
fwhm → smoothing kernel (default: 6 mm)
vox_size → voxel size after normalisation
slice_order → acquisition order (⚠️ confirm with protocol)

Key Features
1. Handles heterogeneous datasets
2. Automatically detects number of slices per session
3. Processes sessions independently (avoids slice mismatch errors)
4. Robust error handling and logging
5. Outputs ready for CONN toolbox analysis

Output
Final preprocessed files have the prefix:
swau*.nii
Where:
u = realigned & unwarped
a = slice timing corrected
w = normalized to MNI
s = smoothed

A log file is saved at:
output_dir/preprocessing_log.txt

Requirements
MATLAB
SPM12 (added to MATLAB path)

Notes
1. Always verify slice order from the MRI acquisition protocol
2. Ensure TR values are correct for each dataset
3. Only 4D NIfTI files (>100MB) are processed as functional data
4. The script processes up to 2 sessions per subject

Use Case
This pipeline is ideal for:
Multi-site studies
Retrospective datasets with inconsistent formats
Large cohort preprocessing for connectivity analysis (e.g., CONN)
