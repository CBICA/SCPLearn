#!/usr/bin/env python
###########################################################################
# @package test_2d_compare
#  @brief This script opens the result of test and checks with ground-truth
# @author Harini Eavani
#
# @Link: https://www.cbica.upenn.edu/sbia/software/
#
# @Contact: sbia-software@uphs.upenn.edu
##########################################################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# first generate 2d scp patterns
# circle
Basis = np.zeros((121,3))

x = np.arange(-5,6,1)
y = np.arange(-5,6,1)
xc,yc = np.meshgrid(x,y)

zc = np.sqrt(xc**2 + yc **2)
zc = np.logical_and(zc > 3.6, zc <4.4 )
zc = zc.astype('int')
Basis[:,0] = zc.ravel()


# square
zc = np.zeros((11,11))
zc[2:8,2] = 1
zc[2:8,8] = 1
zc[2,2:8] = 1
zc[8,2:9] = 1
Basis[:,1] = zc.ravel()

#cross
zc = (xc == yc).astype('int')
zc += (- xc == yc).astype('int')
zc[0,0] = 0
zc[0,10] = 0
zc[10,0] = 0
zc[10,10] = 0
zc = zc.astype('int')
Basis[:,2] = zc.ravel()

norm_factor=np.sqrt(np.sum(Basis*Basis,axis=0))[np.newaxis,:]
Basis=Basis/norm_factor


#load scplearn results
D = 11
T = 100
N = 10
scpFile='test/test_2d_SCPs.mat'
mat_contents = loadmat(scpFile)
B = np.asarray(mat_contents['B'])
norm_factor=np.sqrt(np.sum(B*B,axis=0))[np.newaxis,:]
B=B/norm_factor
innprod = np.zeros(3)
for i in range(3):
    innprod[i] = np.abs(np.dot(B[:,i],Basis[:,i]))

if np.sum(innprod>0.70)==3:
    print('Match between results and ground truth '+str(np.mean(innprod)*100)+'>70% - test passed!')
else:
    print('Match between results and ground truth '+str(np.mean(innprod)*100)+'<70% - test failed!')
    
