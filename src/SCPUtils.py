###########################################################################
# @package SCPUtils
# @brief This script has basic functions for exception handling, execute commands, 
# generate working directories etc \n
#
#
# @author Harini Eavani
#
# <a href="https://www.cbica.upenn.edu/sbia/software/">Link to CBICA Software</a> 
##########################################################################
import sys, subprocess, os
from numpy import  ndarray, floating
from types import ListType

class SbiaException(Exception):
    """our exception class. You can raise it with a msg"""
#  def _get_message(self, message): return self._message
#  def _set_message(self, message): self._message = message
#  message = property(_get_message, _set_message)

    def __init__(self, msg):
        self.message = msg
    def __str__(self):
        return repr(self.message)

def cryandexit(error, filename=None, line=None, exit=True):
    """
    raises an error and exits
      error    - the error text
      filename - the thing the error is associated with
      line     - the line the error was encountered on
      exit     - True if you want this to raise a sys.exit
    """
    if (line == None and filename == None):
        msg='*** Error: %s!' % (error)
    elif (line == None):
        msg='*** Error in %s: %s!' % (str(filename), error)
    else:
        msg='*** Error in %s, line %s: %s!' % (str(filename), str(line), error)
    if exit:
        sys.exit(msg)
    else:
        print msg
        
def execCmd(cmdArray,verbose=0,simulate=False,quiet=False,noFail=False,shell=False):
  """
  executes a command array by space-concatenating it and passing it to subprocess.Popen
  returns a tuple containing (return value, stdout, stderr)
  if the exit status is not 0, raises an SbiaException
    cmdArray - the input command array
    verbose  - turn on extra output
    simulate - simulate the run (no actual commands run)
    quiet    - turns off output of stdout and stderr to the screen
    noFail   - if true does not raise an exception if return value is not 0
    shell    - command can be a bash command
  """
  MAX_CMD_ARR_LEN = 20 # max number of items to dump
  
  cmdArray = [str(i) for i in cmdArray]
  if len(cmdArray) > MAX_CMD_ARR_LEN and verbose<2:
    cmd = " ".join(cmdArray[:MAX_CMD_ARR_LEN]) + " ...[" + str(len(cmdArray)) + "]"
  else:
    cmd = " ".join(cmdArray)
  
  # if verbose, print command which will be run
  if verbose or simulate:
    print '>>>', cmd
  
  # execute command
  ret = [0]
  if simulate:
    ret = [0, 'Simulated', 'Simulated']
  else:
    try:
      # open subprocess
      if shell:
        process = subprocess.Popen(' '.join(cmdArray), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
      else:
        process = subprocess.Popen(cmdArray, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      
      # read and flush stdout until EOF
      while not quiet:
        out = process.stdout.read(1)
        if out == '' and process.poll() != None:
          break
        if out != '':
          sys.stdout.write(out)
          sys.stdout.flush()
    
      # wait until subprocess terminated and set exit code
      err = process.wait()
      # get exit code
      ret[0] = err
      
      # read and flush stderr until EOF
      while not quiet:
        out = process.stderr.read(1)
        if out == '' and process.poll() != None:
          break
        if out != '':
          sys.stderr.write(out)
          sys.stderr.flush()

    except OSError,e:
      msg = cmd + "\tCaused an OSError\n"
      msg = msg + "\t"+e.message
      sys.stderr.write(msg)
      raise e
    except Exception,e:
      sys.stderr.write(e.message)
      raise e

  if verbose>0 or simulate:
    print "return value: %d" % ret[0]

  # if this failed then throw an exception
  if ret[0] != 0 and not noFail:
    msg = cmd + " Failed!"
    raise SbiaException(msg)

  return ret
  
###########################################################################
##########   FILE FUNCTIONS   #############################################
###########################################################################

def fileExists(fileName,isDir=False):
    """
    checks whether a file or dir exists
      filename - path to file
      isDir    - True if you want to check whether it is a directory
    """
    if isDir:
        return(os.path.isdir(fileName))
    else:
        return(os.path.isfile(fileName))  

def isWritable(fileName):
    """
    checks whether directory is writable
    """
    return os.access(fileName,os.W_OK)
        
def getFileBase(fileName,includePath=False,excludeView=False):
    """returns the name part of fileName (the view parameter refers to the Afni +orig suffix)"""
    if includePath:
        root,ext = os.path.splitext(fileName)
    else:
        root,ext = os.path.splitext(os.path.basename(fileName))
    if ext.lower() == '.gz':
        root2,ext2 = os.path.splitext(root)
        if ext2.lower() == '.nii' or ext2.upper() == '.HEAD':
            root = root2
    if excludeView and (root.endswith("+orig") or root.endswith("+acpc") or root.endswith("+tlrc")):
        root = root[:-5]
    return root

def getFileExt(fileName,includeView=False):
    """returns the .extension part of fileName (the view parameter refers to the Afni +orig suffix)"""
    root,ext = os.path.splitext(os.path.basename(fileName))
    if ext.lower() == '.gz':
        root2,ext2 = os.path.splitext(root)
        if ext2.lower() == '.nii' or ext2.upper() == '.HEAD':
            ext = ext2+ext
            root = root2
    if includeView and (root.endswith("+orig") or root.endswith("+acpc") or root.endswith("+tlrc")):
        view = root[-5:]
        ext = view + ext
    return ext

def getFilePath(fileName):
    """returns the path to fileName"""
    return os.path.dirname(fileName)
    

###########################################################################
##########   NUMBER FUNCTIONS   ###########################################
###########################################################################

def isNumber(x, strOk=False):
    """returns True if x is a number (or list of numbers)"""
    if isList(x):
        return len(filter((lambda i : not isNumber(i, strOk=strOk)),x)) == 0
    if strOk: x = tryNumber(x)
    return isinstance(x, (int, long, float, floating))

def isInteger(x):
    """returns True if x is an integer (or list of integers)"""
    if isList(x):
        return len(filter((lambda i : not isInteger(i)),x)) == 0
    return isinstance(x, (int, long))

def isDecimal(x, strOk=False):
    """returns True if x is a float (or list of floats)"""
    if isList(x):
        return len(filter((lambda i : not isDecimal(i)),x)) == 0
    if strOk: x = tryNumber(x)
    return isinstance(x, (float, floating))

def isList(x):
    """returns True if x is some kind of a list"""
    return isinstance(x, (ListType, ndarray, tuple))

def tryNumber(x):
    """tries to return an appropriate number from x or just x if not a number"""
    if isList(x):
        return [tryNumber(i) for i in x]
    if not str(x).isdigit():
        try:
            return float(x)
        except ValueError, e:
            return x
    else:
        try:
            return int(x) # will return long if long int!
        except ValueError, e:
            return x

def extractNumber(x):
    """returns the first integer encountered in a string x (or list of strings) or None if no number is found"""
    if isList(x):
        return [extractNumber(i) for i in x]
    start=None
    x = str(x)
    for i in range(len(x)+1):
        if start!=None and (i==len(x) or not x[i].isdigit()): #return
            return tryNumber(x[start:i])
        elif start==None and x[i].isdigit(): #record the start
            start=i
    return None    
