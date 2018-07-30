#!/usr/bin/env python
###########################################################################
# @package replaceLabels_nib
# @brief This script replaces labels in an atlas with values from a csv column \n
#
# It takes as input \n
# (1) the mask of ROIs/parcellation/atlas \n
# (2) a csv file containing the values that replace the labels \n
# It outputs  \n
# (1) a nifti file with the labels replaced with values
#
#
# @author Harini Eavani
#
# <a href="https://www.cbica.upenn.edu/sbia/software/">Link to CBICA Software</a> 
##########################################################################

import os, sys, getopt
import pandas as pd
from numpy import sort, float32, uint8, zeros, unique
from SCPUtils import *
import nibabel as nib

SVN_REVISION = "1.0.0"
EXEC_NAME = "replaceLabels_nib.py"

def version():
    """prints Version information"""
    msg = EXEC_NAME + """

  Release Information
      replaceLabels_nib.py Release : 1.0.0
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
    Replace labels in an ROI volume with a set of statistical values

  Usage: %(EXEC)s [OPTIONS]

  Required Options:
    [-d --dataFile]   Specify the input spreadsheet (csv format). Must have IDs containing the
                      label number in the first column and descriptive headers in the first row
    [-A --atlas]      Specify the ROI atlas
    [-H --header]     Specify the header of the column whoose values you want to use
                      (optional if your list contains only two columns)
  Other Options:
    [-o --outputDir]    The output directory. Default is the directory of the input spreadsheet
    [-p --prefix]       Specify a prefix for the output file. Default is <atlas prefix>_<header used>

  Get Help:
    [-h --help -u --usage]   Display this message
    [-V --Version]           Version information
    [-v --verbose]           Verbose output

  Examples:
    %(EXEC)s -d roi_stats.csv -A JHU_2m_las_WMPM2.hdr -H GLM_AGE_BETA
   # Uses a transposed spreadsheet, such as that output by statistics.py
    %(EXEC)s -d list.txt -A jakob_rad_convention_labels_no_cere_lps.hdr 
   # Uses a tab-delimited list (created in excel for example)
   # no need for -H if your list contains only 2 columns
  """ % {'EXEC':EXEC_NAME}

def main():
    """parse command line input and call appropriate functions"""
    verbose  = 0

    dataFile   = None
    atlas      = None

    prefix     = None
    outDir     = None
    header     = None

    rOpts = 0
    try:
        opts, files = getopt.gnu_getopt(sys.argv[1:], "hd:o:p:vVuA:H:",
        ["help","dataFile=","outputDir=","prefix=","verbose","Version","usage","atlas=","header="])

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
        elif o in ["-H", "--header"]:
            header = a
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

    # grab the input spreadsheet
    mydata = pd.read_csv(dataFile,header=0,quotechar='"',sep=',',na_values = ['NA', '-', '.', ''],index_col=0) 
    # check size
    if mydata.shape[1] == 0:
        cryandexit("No columns found in your spreadsheet! Wrong delimiter?")
    elif mydata.shape[1] < 1:
        cryandexit("Too few columns found in your spreadsheet! Need at least 2")
    if mydata.shape[0] == 0:
        cryandexit("No rows found in your spreadsheet!")
    elif mydata.shape[0] < 1:
        cryandexit("Too few rows found in your spreadsheet! Need at least 2")

    # grab the only header available
    if not header:
        if mydata.shape[0] == 1:
            header = mydata.keys()[0]
        else:
            cryandexit("More than two columns found! Please specify -H")

    # build output filename
    if not prefix:
        prefix = getFileBase(atlas) + '_' + header.replace('__','')
    outputfile = os.path.join(outDir,prefix+getFileExt(atlas))

    # build dictionary
    mydict = {}
        # check if header is there
    if not header in mydata.columns:
        cryandexit("Unable to find column in spreadsheet",header)
    myids = mydata.index.values 
    # build replacement dictionary
    for id in myids:
        replacement = mydata.loc[id][header]
        if not replacement:
            replacement = 0
        elif isNumber(replacement,strOk=True):
            replacement = tryNumber(replacement)
        else:
            raise SbiaException("Illegal value encountered at " + str(id) + ", " + str(header))
        #mydict[extractNumber(id+1)] = replacement
        mydict[id] = replacement
    
    
    atlas_im = nib.load(atlas)
    atlas_im.set_data_dtype(float32)
        
    atlas_data = atlas_im.get_data()
    out_data = zeros(atlas_data.shape)
    labels = sort(unique(atlas_data.flatten())) # note labels begin from 1

    # replace labels!
    for label in labels:
        if not mydict.has_key(label):
            if verbose:
                print "!! Warning: Could not find replacement for label",label,".Setting it to 0"
            out_data[atlas_data==label] = 0
        else:
            out_data[atlas_data==label] = mydict[label]

    # save
    out_img = nib.Nifti1Image(out_data, atlas_im._affine,header=atlas_im._header)
    nib.save(out_img,outputfile)
    if verbose:
        print ">> Saved",outputfile

    sys.exit(0)


if __name__ == '__main__': main()