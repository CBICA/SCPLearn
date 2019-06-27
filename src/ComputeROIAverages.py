#!/usr/bin/env python
###########################################################################
# @package ComputeROIAverages
# @brief This script extracts average ROI time-series from nifti data \n
#
# It takes as input \n
# (1) list of nifti files \n
# (2) the mask of ROIs/parcellation/atlas \n
# It outputs  \n
# (1) mat files containing all the outputs in ROI space  \n
# (2) a text file containing the list of all mat files that were generated
#
#
# @author Harini Eavani
#
# <a href="https://www.cbica.upenn.edu/sbia/software/">Link to CBICA Software</a> 
##########################################################################

import os, sys
from numpy import sort, uint8, unique, mean
from SCPUtils import *
import nibabel as nib
#import time
from scipy.io import savemat


   
def ComputeROIAverages(dataFile,atlas,outDir):
    
    dataFile = os.path.realpath(dataFile)
    if not fileExists(dataFile):
        cryandexit("File does not exist", dataFile)
    
    atlas = os.path.realpath(atlas)
    if not fileExists(atlas):
        cryandexit("File does not exist", atlas)
            
    
    outDir = os.path.realpath(outDir)
    if not os.path.exists(outDir):
        print "!! Output path does not exist. Creating", outDir
        makeDir(outDir)
    if not isWritable(outDir):
        cryandexit("Output directory is not writable", outDir)
 

    atlas_im = nib.load(atlas)
    atlas_im.set_data_dtype(uint8)    
    atlas_data = atlas_im.get_data()

    print('Loaded atlas data')
    labels = sort(unique(atlas_data.flatten()))
    labels = list(set(labels) - set([0]))
    
    print('Identified unique labels')
    L = len(labels)
    
    # label indexes
    label_indexes = [atlas_data==labels[e] for e in xrange(L)]

    # load subject data
    list_of_files = os.path.join(outDir,getFileBase(dataFile) + '_' + getFileBase(atlas)+'.txt')
    with open(list_of_files, "wb") as flist:
        with open(dataFile, "r") as fin:
            for subjectFile in fin:
                prefix = getFileBase(subjectFile) + '_' + getFileBase(atlas)
                outputfile = os.path.join(outDir,prefix+'.mat')
                if not fileExists(outputfile):
                    print("Roi time-series for  "+ subjectFile+ " not found, generating..")
                    subject_im = nib.load(subjectFile.rstrip('\n'))
                    subject_data = subject_im.get_data()
                    print('Loaded subject '+subjectFile)
                
                    # check if both volumes have same dimensions
                    if not (atlas_data.shape[0:3] == subject_data.shape[0:3]):
                        cryandexit("Subject and Atlas images have different voxel dimensions. Are they in same space?")
                        
                    #start = time.time()
                    ts = [mean(subject_data[label_indexes[x],:],axis=0) for x in xrange(L)]
                    #end = time.time()
                    #print end - start
                    
                    savemat(outputfile, mdict={'ts': ts},appendmat=False)
                flist.write(outputfile+'\n')
                
    flist.close()
    fin.close()
