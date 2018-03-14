#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 22:04:50 2018

@author: zhi
"""


import pickle
import numpy as np
import json
import time
import pandas

def localMaximumQCM(cMatrix, gamma, t, lambda_val):
    C = []
    # Step 1 - find the local maximal edges
    maxV = cMatrix.max(1)
    maxInd = np.argmax(cMatrix, axis=1)
    maxEdges = []
    maxW = []
    for i in range(0,len(maxInd)):
        if maxV[i] == max(cMatrix[maxInd[i], :]):
            maxEdges.append([maxInd[i], i])
            maxW.append(maxV[i])
    sortMaxV = np.sort(maxW)[::-1]
    sortMaxInd = np.argsort(maxW)[::-1] # orders may different from Matlab results
    sortMaxEdges = np.asarray(maxEdges)[sortMaxInd,:]
    
    currentInit, noNewInit = 0, False
    nodesInCluster = []
    while currentInit <= len(sortMaxInd) - 1 and noNewInit == False:
        if sortMaxV[currentInit] < gamma * sortMaxV[0]:
            noNewInit = True
        else:
            if (sortMaxEdges[currentInit, 0] in nodesInCluster) == 0 \
                    and (sortMaxEdges[currentInit, 1] in nodesInCluster) == 0:
                newCluster = sortMaxEdges[currentInit, :]
                addingMode = 1
                currentDensity = sortMaxV[currentInit]
                nCp = 2
                totalInd = np.asarray(range(0, cMatrix.shape[0]))
                remainInd = np.setdiff1d(totalInd, newCluster)
                while addingMode == 1:
                    neighborWeights = np.sum(cMatrix[np.ix_(newCluster, remainInd)], axis=0)
                    maxNeighborWeight = neighborWeights.max()
                    maxNeighborInd = np.argmax(neighborWeights)
                    c_v = maxNeighborWeight/nCp
                    alphaN = 1 - 1./(2*lambda_val*(nCp+t))
                    if (c_v >= alphaN * currentDensity):
                        newCluster = np.append(newCluster, remainInd[maxNeighborInd])
                        nCp += 1
                        currentDensity = (currentDensity*((nCp-1)*(nCp-2)/2.)+maxNeighborWeight)/(nCp*(nCp-1)/2.)
                        remainInd = np.setdiff1d(remainInd, remainInd[maxNeighborInd])
                    else:
                        addingMode = 0
                nodesInCluster = np.concatenate((nodesInCluster, newCluster)).astype(int)
                C.append(newCluster.tolist())
        currentInit += 1
    return C

def lmQCMMerge(C, beta, minClusterSize):
    sortC = sorted(C,key=len)[::-1]
    mergedCluster = filter(lambda x: len(x) >= minClusterSize, sortC)
    mergeOccur, currentInd = 1, -1
    print "start merge"
    while mergeOccur == 1:
        print "merging ..."
        mergeOccur = 0
        while currentInd < len(mergedCluster):
            currentInd += 1
            excludeInd = []
            if currentInd < len(mergedCluster):
                keepInd = range(0, currentInd+1)
                for j in range(currentInd+1, len(mergedCluster)):
                    j = 1
                    
                    interCluster = np.intersect1d(mergedCluster[currentInd], mergedCluster[j])
                    
                    # Known issue:
                    # Sorting gives different mergedCluster in python and Matlab.
                    # Due to the ordering of the mergedCluster, say mergedCluster is with size [{1x5},{1x4},{1x4},{1x3}]
                    # but in matlab, the 2nd appeared {1x4} cluster may be the 1st appeared {1x4} cluster, which will
                    # cause different result in this section.
                    print(j)
                    print(beta*min(len(mergedCluster[j]), len(mergedCluster[currentInd])))
                    print(len(interCluster))
                    print(len(mergedCluster[currentInd]))
                    
                    if len(interCluster) >= beta*min(len(mergedCluster[j]), len(mergedCluster[currentInd])):
                        mergedCluster[currentInd] = np.union1d(mergedCluster[currentInd], mergedCluster[j]).tolist()
                        mergeOccur = 1
                    else:
                        keepInd.append(j)
                    j = j+1
                mergedCluster = [ mergedCluster[i] for i in keepInd ]
                print "merged cluster length: " + str(len(mergedCluster))
        sizeMergedCluster = np.zeros(len(mergedCluster))
        for i in range(0,len(mergedCluster)):
            sizeMergedCluster[i] = len(mergedCluster[i])
        sortSize = np.sort(sizeMergedCluster)[::-1]
        sortMergedInd = np.argsort(sizeMergedCluster)[::-1]
        mergedCluster = [ mergedCluster[i] for i in sortMergedInd ]
        currentInd = 0
    return mergedCluster







with open('/Users/zhi/Desktop/inputdata.pkl', 'rb') as f:
   inputdata = pickle.load(f)
nRow, nCol = 2000, 415
gamma,t,lambda_val,beta,minClusterSize = 0.55,1,1,0.4,10
finalExpArray = inputdata
method = "Pearson"
finalExp = np.asarray(finalExpArray).reshape((nCol, nRow))
df = pandas.DataFrame(finalExp)

if method == "Pearson":
	cMatrix = np.asarray(df.corr("pearson"))
if method == "Spearman":
	cMatrix = np.asarray(df.corr("spearman"))
    
np.fill_diagonal(cMatrix, 0)
    
# starting lmQCM
tm = time.time()
cMatrix = abs(cMatrix)
cell = localMaximumQCM(abs(cMatrix), gamma, t, lambda_val)
print "step 2: calculate lmQCM elapsed time: " + str("{0:.2f}".format(time.time() - tm)) + " seconds"

# step 3: Merge the overlapped networks
tm = time.time()
mergedCluster = lmQCMMerge(cell, beta, minClusterSize)
print "step 3: Merge the overlapped networks elapsed time: " + str("{0:.2f}".format(time.time() - tm)) + " seconds"
print "Length of merged cluster: " + str(len(mergedCluster))






