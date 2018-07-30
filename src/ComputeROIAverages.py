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

import os, sys, getopt
import csv
from numpy import sort, uint8, unique, mean
from SCPUtils import *
import nibabel as nib
import time
from scipy.io import savemat

SVN_REVISION = "1.0.0"
EXEC_NAME = "replaceLabels.py"

def version():
    """prints Version information"""
    msg = EXEC_NAME + """

  Release Information
      ComputeROIAverages.py Release : 1.0.0
      Contact : SBIA Group <sbia-software at uphs.upenn.edu>
    """

    svnMsg = """
  SVN information
      Project Revision : """ + SVN_REVISION

    if (len(SVN_REVISION) != 0):
        msg = msg + svnMsg
    print msg


def usage():
    """prints usage information"""
    print r"""
  %(EXEC)s--
    Compute average ROI values given an atlas.

  Usage: %(EXEC)s [OPTIONS]

  Required Options:
    [-d --dataFile]   Specify the list of input NIFTI images. 
    [-A --atlas]      Specify the ROI atlas, also as NIFTI image.

  Other Options:
    [-o --outputDir]    The output directory. Default is the directory of the input spreadsheet

  Get Help:
    [-h --help -u --usage]   Display this message
    [-V --Version]           Version information
    [-v --verbose]           Verbose output

  Examples:
    %(EXEC)s -d projname_all_subjects.txt -A JHU_2m_las_WMPM2.hdr 
  """ % {'EXEC':EXEC_NAME}
  
   
def main():
    """parse command line input and call appropriate functions"""
    verbose  = 0

    dataFile   = None
    atlas      = None

    prefix     = None
    outDir     = None


    rOpts = 0
    try:
        opts, files = getopt.gnu_getopt(sys.argv[1:], "hd:o:p:vVuA:",
        ["help","dataFile=","outputDir=","prefix=","verbose","Version","usage","atlas="])

    except getopt.GetoptError, err:
        usage()
        cryandexit(str(err)) # option -x not recognized..

    for o, a in opts:
        #print o,a
        if o in ["-v", "--verbose"]:
            verbose+=1
        elif o in ["-h", "--help","-u","--usage"]:
            usage()
            sys.exit(0)
        elif o in ["-V", "--Version"]:
            version()
            sys.exit(0)
        elif o in ["-d", "--dataFile"]:
            dataFile = os.path.realpath(a)
            if not fileExists(dataFile):
                cryandexit("File does not exist", dataFile)
            rOpts+=1
        elif o in ["-p", "--prefix"]:
            prefix = a
        elif o in ["-o", "--outputDir"]:
            outDir = a
        elif o in ["-A", "--atlas"]:
            atlas = os.path.realpath(a)
            if not fileExists(atlas):
                cryandexit("File does not exist", atlas)
            rOpts+=1
        else:
            assert False, "unhandled option"

    if rOpts != 2:
        usage()
        cryandexit("Please specify all required options")

    # parameter checking
    if not outDir:
        outDir = getFilePath(dataFile)
    outDir = os.path.realpath(outDir)
    if not os.path.exists(outDir):
        print "!! Output path does not exist. Creating", outDir
        makeDir(outDir)
    if not isWritable(outDir):
        cryandexit("Output directory is not writable", outDir)
 

    atlas_im = nib.load(atlas)
    atlas_im.set_data_dtype(uint8)    
    atlas_data = atlas_im.get_data()
    if verbose:
        print('Loaded atlas data')
    labels = sort(unique(atlas_data.flatten()))
    labels = list(set(labels) - set([0]))
    
    if verbose:
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
                    print("Roi time-series for  ", subjectFile," not found, generating..")
                    subject_im = nib.load(subjectFile.rstrip('\n'))
                    subject_data = subject_im.get_data()
                    if verbose:
                        print('Loaded subject '+subjectFile)
                
                    # check if both volumes have same dimensions
                    if not (atlas_data.shape[0:3] == subject_data.shape[0:3]):
                        cryandexit("Subject and Atlas images have different voxel dimensions. Are they in same space?")
                        
                    start = time.time()
                    ts = [mean(subject_data[label_indexes[x],:],axis=0) for x in xrange(L)]
                    end = time.time()
                    print end - start
                    
                    savemat(outputfile, mdict={'ts': ts},appendmat=False)
                flist.write(outputfile+'\n')
                
    flist.close()
    fin.close()
    
    if verbose:
        print ">> Saved",outputfile

    sys.exit(0)


if __name__ == '__main__': main()