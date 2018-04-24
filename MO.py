# -*- coding: utf-8 -*-
"""
Spyder Editor

Script reads Gaussian MO coefficients from a file as:
1a1b2a2b
and outputs them as
1a2a1b2b
"""
import numpy as np
import math

def deInterlace(inMatrix, dimen):
    tempMatrix1 = np.ndarray(shape=(dimen,dimen))
    tempMatrix2 = np.ndarray(shape=(dimen,dimen))
    nbasis = int(dimen/2)
    
    "Deinterlace rows"
    for i in range(0,dimen):
        rowoff = 0
        for j in range(0,nbasis):
            tempMatrix1[j,i] = inMatrix[j+rowoff,i]
            tempMatrix1[j+nbasis,i] = inMatrix[j+1+rowoff,i]
            rowoff += 1
        
    
    "Deinterlace columns"
    coloff = 0
    for i in range(0,nbasis):
        for j in range(0,dimen):
            tempMatrix2[j,i] = tempMatrix1[j,i+coloff]
            tempMatrix2[j,i+nbasis] = tempMatrix1[j,i+1+coloff]
        coloff += 1
        
    return tempMatrix2

def deInterlaceRows(inMatrix, dimen):
    tempMatrix1 = np.ndarray(shape=(dimen,dimen))
    nbasis = int(dimen/2)
    
    "Deinterlace rows"
    for i in range(0,dimen):
        rowoff = 0
        for j in range(0,nbasis):
            tempMatrix1[j,i] = inMatrix[j+rowoff,i]
            tempMatrix1[j+nbasis,i] = inMatrix[j+1+rowoff,i]
            rowoff += 1
            
    return tempMatrix1
    

"Script main begins"
MOlist = []
f = open('O2_triplet','r')
for line in f:
    MOlist += line.split()
    
elements = len(MOlist)
dimension = int(math.sqrt(elements/2))
nbasis = int(dimension/2)
print(nbasis)

realCoeff = np.ndarray(shape=(dimension,dimension))
imagCoeff = np.ndarray(shape=(dimension,dimension))

for i in range(dimension):
    for j in range(dimension):
        realCoeff[i,j] = 0.0
        

"propogate coeff matricies:"
iteration = 0
for i in range(0,elements-1,int(dimension*2)): 
    offset = 0
    for j in range(0,dimension):
        realCoeff[j,iteration] = MOlist[i+offset]
        imagCoeff[j,iteration] = MOlist[i+1+offset]
        offset +=2
        
    iteration += 1

realMQC = deInterlaceRows(realCoeff,dimension)
imagMQC = deInterlaceRows(imagCoeff,dimension)

    
for i in range(dimension):
    stringout = ''
    for j in range(dimension):
        stringout += '  {0: 5f}'.format(realCoeff[i,j])
    print(stringout)
    
print("\n")

for i in range(dimension):
    stringout = ''
    for j in range(dimension):
        stringout += '  {0: 5f}'.format(imagCoeff[i,j])
    print(stringout)
    
print("\n")
    
for i in range(dimension):
    stringout = ''
    for j in range(dimension):
        stringout += '  {0: 5f}'.format(realMQC[i,j])
    print(stringout)
    
print("\n")
    
for i in range(dimension):
    stringout = ''
    for j in range(dimension):
        stringout += '  {0: 5f}'.format(imagMQC[i,j])
    print(stringout)
    


    
    
    

    
