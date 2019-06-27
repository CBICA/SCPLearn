\mainpage SCPLearn 

  Section of Biomedical Image Analysis <br>
  Department of Radiology <br>
  University of Pennsylvania <br>
  Goddard Building <br>
  3701 Hamilton Walk, 6th Floor <br>
  Philadelphia, PA 19104 <br>

  Web:   http://www.cbica.upenn.edu/sbia/ <br>
  Email: sbia-software at uphs.upenn.edu

  Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
  See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.

<b> Author: </b>
Harini Eavani
Harini.Eavani@uphs.upenn.edu 
cbica-software@uphs.upenn.edu

===============
1. INTRODUCTION
===============

This software generates Sparse Connectivity Patterns from resting state fmri connectivity data.
It applies sparse matrix decomposition to rsfMRI connectivity matrices to obtain group-common 
Sparse Connectivity Patterns (SCPs) and subject-specific SCP coefficients. SCP coefficients encode 
the total connectivity of a pattern within a subject. Each SCP can be further sub-divided into smaller SCPs,
along with their associated subject level coefficients. <br>

For details see the following paper:
http://www.sciencedirect.com/science/article/pii/S1053811914008003

===============
2. SOFTWARE
===============

This software has been tested on Linux operating systems only. It is implemeted primarily in MATLAB. For command-line argument parsing, and nifti I/O operations, python code is used.

==========
3. INSTALL    
==========

----------------
3.1 Dependencies
----------------
- MATLAB Compiler mcc version 5.2 (R2014B)
- MATLAB R2014B
- Python 2.7.9
- Python library numpy 1.7.2
- Python library scipy 0.15.1
- Python library pandas 0.16.2
- Python library nibabel 2.0.1

Make sure all dependencies are met before proceeding with install.

-------------------------------------
3.2 Generating standalone executables
-------------------------------------
1)  Within the <code>src/</code> directory, run make:

<code>./make</code>

2)  Within the <code>src/</code> directory, test the code as follows:

<code>./make test</code>

This test runs scplearn on sythetic test images, and compares the results to the ground-truth.
This command should return three SCPs, saved as nifti files in <code> src/test/test_2d_SCP_{1,2,3}.nii.gz </code> along with other output. <br>

By design, three patterns were used as 'ground-truth' SCPs to generate the synthetic data; a cross, a circle and a square. <br>
Was SCPLearn able to clearly separate these three patterns in the output? You can compare your result to the results in <code> src/test/Test_2d_all.png </code>

If testing succeeded, you will see the message: 
"Match between results and ground truth xx.xxxx>70% - test passed!", where xx.xx is the percentage match between result and ground-truth.
Proceed to install only if testing is completed.

3) Within the <code>src/</code> directory, run make install with the location of the install directory as the argument:

<code>./make install <installDir></code>

4) Add the install directory to your path by running the following command. 

Replace <code>${installDir}</code> with the location of the install directory from step 1 above.

<code>export PATH=${PATH}:${installDir}</code>

----------------------------
3.3 Generating code documentation
----------------------------

The above make install command creates doxygen documentation in "docs/*". Please open "docs/index.html" in you browser.

==========
4. Usage
==========

This software provides two executibles for users - scplearn and scptest.

scplearn takes a set of data as input(either NIFTI data or mat files), and generates both SCPs as well as SCP coeffcients for that data. SCPs are common to the whole group, SCP coefficients are subject specific.

	scplearn--
    Generates Sparse Connectivity Patterns from resting state fMRI data

    Usage: scplearn [OPTIONS] 

    Required Options:
    [-d --data]    Specify the text file with list of nifti files (required) 
                    *** OR ***
                   Specify the text file with list of mat files (required) 
                   Each mat file must contain a variable named 'ts'
                   'ts' must be a matrix of time-series, size (# of ROIs X # of timepoints)
                    
    [-m --mask]    Specify the nifti ROI/parcel/atlas file (required)
    [-p --prefix]      Specify the prefix of the output file  (required)
    
    Options:
    [-t --type]   Specify the data type of input with either "matlab" or "nii".      
    [-s --sparsity]   Specify the sparsity constraint as positive value. Default = nROIs/10. 
    [-n --numberOfSCPs]   Specify number of primary SCPs. Default = 10.    
    [-r --pruning]   Specify the pruning threshold as a value between [0,1]. Default = 0.7.    
    [-l --levels]   Specify number of secondary SCPs. Default = 50.

    [-o --outputDir]       The output directory to write the results. Defaults to the location of the input file
    [-w --workingDir]      Specify a working directory. By default a tmp dir is created and used
    [-u --usage | -h --help]    Display this message
    [-v --verbose]              Verbose output
    [-V --Version]              Display version information

    Examples:
    scplearn -d list_of_nifti_files.txt -t nii -m Grasp_level5.nii -p SCP_results -n 10 -s 50 -o /sbia/sbiaprj/BLSA -v
    scplearn -d list_of_mat_files.txt -m Grasp_level5.nii -p SCP_results -n 10 -r 0.5 -o /sbia/sbiaprj/BLSA -v
    scplearn -d list_of_mat_files.txt -t matlab -m Grasp_level5.nii -p SCP_results -l -n 10 -o /sbia/sbiaprj/BLSA -v
    
    Example list_of_nifti_files.txt:
    ProjName_subj_165464.nii.gz
    ProjName_subj_26464.nii.gz
    ProjName_subj_1054.nii.gz......      

    

scplearn takes as input a text file with the list of nifti/mat files. If input is nifti, each nifti file must be pre-processed rsfMRI 4D 
scan in a common template/standard space. If input is mat files, each mat file must contain a variable named 'ts'. 'ts' must be a matrix of time-series, 
size (# of ROIs X # of timepoints) It also requires the "node" definitions (e.g. atlas, parcellation, ROI definition, etc) to be input as a nifti file. 
The node definitions need to be in the same template space as the other nifti files.

Mandatory arguments are:
 - List of nifti/mat files
 - NIFTI mask/atlas/ROI file
 - Prefix for all output files

Optional arguments are:
 - Number of SCPs at the primary level. Default = 10.
 - Sparsity level of SCPs. Default = Number of ROIs / 10. 
 - Pruning threshold. Default = 0.7. SCPs with inner-product overlap greater than this threshold are discarded
 - Number of SCPs at the secondary level in the hierarchical decomposition. Default=50. 
 
In the output directory, the software returns:
 - One NIFTI file for each SCP that is generated, with file name <code><prefix>_SCP_#.nii.gz</code>
 - A <code><prefix>_SCP_Coeffs.csv</code> file with the SCP coefficients for all the subjects, indexed by the NIFTI/mat filename that was input
 - A <code><prefix>_SCP_basis.csv</code> file with the SCP basis, in csv format
 - A <code><prefix>_SCPs.mat</code> file with all the outputs in MATLAB '.mat' format
 - A <code><prefix>_ts.mat</code> file with the average time-series from all the subjects in MATLAB '.mat' format
 - A <code><subject>_<atlas>.mat</code> file with the average time-series for each subject 


scptest takes a set of data AND SCPs as input, and generates on SCP coefficients for the data, corresponding to the SCPs that were input.
For example, scplearn may be used when generating SCPs for the first time within a study. For subsequent participants that are acquired and processed, there is no need to re-compute SCPs; one can use existing SCPs as the reference and generate only SCP coefficients for the subsequent participants.
 
    scptest--
    Computes SCP coefficients for input data given input SCP basis

    Usage: scptest [OPTIONS] 

    Required Options:
    [-d --data]    Specify the text file with list of nifti files (required) 
                    *** OR ***
                   Specify the text file with list of mat files (required) 
                   Each mat file must contain a variable named 'ts'
                   'ts' must be a matrix of time-series, size (# of ROIs X # of timepoints)
    [-m --mask]    Specify the nifti ROI/parcel/atlas file (required)                
    [-b --basis]    Specify the input SCP basis, saved as variable 'B' in mat file (required)
    [-p --prefix]      Specify the prefix of the output file  (required)
    
    Options:
    [-o --outputDir]       The output directory to write the results. Defaults to the location of the input file
    [-w --workingDir]      Specify a working directory. By default a tmp dir is created and used
    [-u --usage | -h --help]    Display this message
    [-v --verbose]              Verbose output
    [-V --Version]              Display version information

    Examples:
    scptest -d list_of_nifti_files.txt -b BLSA_SCP_Basis.mat -m Grasp_level5.nii -p SCP_test_coeffs -o /sbia/sbiaprj/BLSA -v
    scptest -d list_of_mat_files.txt -b BLSA_SCP_Basis.mat -m Grasp_level5.nii -p SCP_test_coeffs -o /sbia/sbiaprj/BLSA -v
    
    Example list_of_nifti_files.txt:
    ProjName_subj_165464.nii.gz
    ProjName_subj_26464.nii.gz
    ProjName_subj_1054.nii.gz......

scptest takes as input a text file with the list of nifti/mat files, NIFTI "node" file as well as the SCPs stored in "mat" file format.

Mandatory arguments are:
 - List of nifti/mat files
 - NIFTI mask/atlas/ROI file
 - SCPs saved in mat file format (generated by a prior run of scplearn)
 - Prefix for all output files

In the output directory, the software returns:
 - A <code><prefix>_SCP_Coeffs.csv</code> file with the SCP coefficients for all the subjects, indexed by the NIFTI/mat filename that was input
 - A <code><prefix>_SCPs.mat</code> file with all the outputs in MATLAB '.mat' format
 - A <code><prefix>_ts.mat</code> file with the average time-series from all the subjects in MATLAB '.mat' format
 - A <code><subject>_<atlas>.mat</code> file with the average time-series for each subject 

 
===========
6. LICENSING
===========

  See http://www.cbica.upenn.edu/sbia/software/license.html or "licences/COPYING" file.
