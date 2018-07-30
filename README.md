\mainpage SCPLearn 

  Section of Biomedical Image Analysis <br>
  Department of Radiology <br>
  University of Pennsylvania <br>
  Goddard Building <br>
  3701 Hamilton Walk, 6th Floor <br>
  Philadelphia, PA 19104 <br>

  Web:   http://www.cbica.upenn.edu/sbia/ <br>
  Email: sbia-software at uphs.upenn.edu

  Copyright (c) 2015 University of Pennsylvania. All rights reserved. <br>
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
1) Within the <code>src/</code>  directory, run <code>makeAll.py</code>  with the location of the install directory as the argument:

<code>makeAll.py <installDir></code> 

2) Add the install directory to your path by running the following command. 

Replace <code>${installDir}</code> with the location of the install directory from step 1 above.

<code>export PATH=${PATH}:${installDir}</code> 

----------------------------
3.3 Generating documentation
----------------------------

Run <code>"doxygen Doxyfile"</code> from the parent directory.
The documentation is built in <code>"docs/*"</code>. Open <code>"docs/html/index.html"</code> in your favorite browser to begin.

==========
4. TEST    
==========

To test the software, we have provided synthetic nifti as well as mat data for ten 'subjects'. These are located in <code>src/test</code>. <br>

To run the software on the synthetic data, run any of the following the following command from the <code>src/</code> folder <br>

<code> SCPLearn.py -d test/nifti_list.txt -p test_2d -m test/test_2d_mask.nii.gz -o test/ -n 3 -s 10 </code> <br>

or

<code> SCPLearn.py -d test/mat_list.txt -p test_2d -m test/test_2d_mask.nii.gz -o test/ -n 3 -s 10 </code> <br>


This command should return three SCPs, saved as nifti files in <code> src/test/test_2d_SCP_{1,2,3}.nii.gz </code> along with other output. <br>

By design, three patterns were used as 'ground-truth' SCPs to generate the synthetic data; a cross, a circle and a square. <br>
Was SCPLearn able to clearly separate these three patterns in the output? You can compare your result to the results in <code> src/test/Test_2d_all.png </code>


==========
5. Usage
==========

    Usage: SCPLearn.py [OPTIONS] 

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
    [-s --sparsity]   Specify the sparsity constraint as positive value. Default = nROIs/8 
    [-n --numberOfSCPs]   Specify the number of SCPs as a number. Default = 10.    
    [-r --pruning]   Specify the pruning threshold as a value between [0,1]. Default = 0.7.    
    [-l --levels]   Hierarchical learning of SCPs. Not run by default. 

    [-o --outputDir]       The output directory to write the results. Defaults to the location of the input file
    [-w --workingDir]      Specify a working directory. By default a tmp dir is created and used
    [-u --usage | -h --help]    Display this message
    [-v --verbose]              Verbose output
    [-V --Version]              Display version information

    Example:
    SCPLearn.py -d list_of_nifti_files.txt -t nii -m Grasp_level5.nii -p SCP_results -n 10 -s 50 -o /sbia/sbiaprj/BLSA -v
    SCPLearn.py -d list_of_mat_files.txt -m Grasp_level5.nii -p SCP_results -n 10 -r 0.5 -o /sbia/sbiaprj/BLSA -v
    SCPLearn.py -d list_of_mat_files.txt -t matlab -m Grasp_level5.nii -p SCP_results -l -n 10 -o /sbia/sbiaprj/BLSA -v

    Example list_of_nifti_files.txt without sample weights:
    ProjName_subj_165464.nii.gz
    ProjName_subj_26464.nii.gz
    ProjName_subj_1054.nii.gz......

    Example list_of_nifti_files.txt with sample weights:
    ProjName_subj_165464.nii.gz,172
    ProjName_subj_26464.nii.gz,160
    ProjName_subj_1054.nii.gz,120......           
    

SCPLearn software takes as input a text file with the list of nifti/mat files. If input is nifti, each nifti file must be pre-processed rsfMRI 4D 
scan in a common template/standard space. If input is mat files, each mat file must contain a variable named 'ts'. 'ts' must be a matrix of time-series, 
size (# of ROIs X # of timepoints) It also requires the "node" definitions (e.g. atlas, parcellation, ROI definition, etc) to be input as a nifti file. 
The node definitions need to be in the same template space as the other nifti files.

Optional arguments are:
 - Number of SCPs at the primary level. Default = 10.
 - Sparsity level of SCPs. Default = Number of ROIs / 8. 
 - Pruning threshold. Default = 0.7. SCPs with inner-product overlap greater than this threshold are discarded
 - Hierarchical decomposition? Default: Not used. 
 
In the output directory, the software returns:
 - One NIFTI file for each SCP that is generated, with file name <code><prefix>_SCP_#.nii.gz</code>
 - A <code><prefix>_SCP_Coeffs.csv</code> file with the SCP coefficients for all the subjects, indexed by the NIFTI filename that was input
 - A <code><prefix>_SCP_basis.csv</code> file with the SCP basis, in csv format
 - A <code><prefix>_SCPs.mat</code> file with all the outputs in MATLAB '.mat' format
 - A <code><prefix>_ts.mat</code> file with the average time-series from all the subjects in MATLAB '.mat' format
 - A <code><subject>_<atlas>.mat</code> file with the average time-series for each subject 

 
 
===========
6. LICENSING
===========

  See http://www.cbica.upenn.edu/sbia/software/license.html or "licences/COPYING" file.
