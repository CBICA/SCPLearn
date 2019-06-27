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

import os
import pandas as pd
from numpy import sort, float32, zeros, unique
from SCPUtils import *
import nibabel as nib

def replaceLabels_nib(dataFile,atlas,header,prefix,outDir):
    
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
            print "!! Warning: Could not find replacement for label",label,".Setting it to 0"
            out_data[atlas_data==label] = 0
        else:
            out_data[atlas_data==label] = mydict[label]

    # save
    out_img = nib.Nifti1Image(out_data, atlas_im._affine,header=atlas_im._header)
    nib.save(out_img,outputfile)
