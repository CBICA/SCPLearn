#!/usr/bin/env python
###########################################################################
# @package makeAll
#  @brief This script calls mcc within matlab 
# It takes as input 
# (1) install directory
# It outputs 
# (1) mcc executible SCPLearnFromMatFiles
# (2) copies all py files to install directory
#
# @Dependencies
#
# @author Harini Eavani
#
# @Link: https://www.cbica.upenn.edu/sbia/software/
#
# @Contact: sbia-software@uphs.upenn.edu
##########################################################################

import os, sys
from SCPUtils import *

if len(sys.argv) < 2:
    cryandexit("Please specify a install directory")

installDir = sys.argv[1]
if not os.path.exists(installDir):
    os.makedirs(installDir)
elif not os.path.isdir(installDir):
    cryandexit("Install dir is not a directory", installDir)

#matlab -nodesktop -nosplash -r "mcc -mv SCPLearnFromMatFiles.m -R -singleCompThread -R -nojvm -d ${installDir}; exit"
cmdArray=['matlab','-nodesktop','-nosplash','-r','"','mcc','-mv','SCPLearnFromMatFiles.m','-R',
              '-singleCompThread','-R','-nojvm','-d',installDir,';','exit','"']
             
execCmd(cmdArray,verbose=1,simulate=False,shell=True)
exeFile=installDir+'/'+'SCPLearnFromMatFiles'
if not fileExists(exeFile):
    cryandexit("Error in building mcc executible, ", exeFile," not found")

#cp *py ${installDir}/
for pyf in ['SCPLearn.py','ComputeROIAverages.py','replaceLabels_nib.py','SCPUtils.py']:
    cmdArray=['cp',pyf,installDir]
    execCmd(cmdArray,verbose=1,simulate=False,shell=True)
    pyFile=installDir+'/'+pyf
    if not fileExists(pyFile):
        cryandexit("Error in copying python files to install dir, ", pyFile," not found")
    #chmod u+x ${installDir}/*
    cmdArray=['chmod','u+x',installDir+'/'+pyf]
    execCmd(cmdArray,verbose=1,simulate=False,shell=True)

# exporting path
# For linux users run
# export PATH=${PATH}:${installDir}
##

