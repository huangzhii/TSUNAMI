#!/Users/zhi/anaconda2/bin/python

import numpy as np
import json
import time
import pandas

def massivePCC(finalExpArray, nRow, nCol, method):
	finalExp = np.asarray(finalExpArray).reshape((nCol, nRow))
	df = pandas.DataFrame(finalExp)
	if nRow > 1:
		if method == "Pearson":
			return np.asarray(df.corr("pearson"))
		if method == "Spearman":
			return np.asarray(df.corr("spearman"))

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
                    interCluster = np.intersect1d(mergedCluster[currentInd], mergedCluster[j])
                    if len(interCluster) >= beta*min(len(mergedCluster[j]), len(mergedCluster[currentInd])):
                        mergedCluster[currentInd] = np.union1d(mergedCluster[currentInd], mergedCluster[j]).tolist()
                        mergeOccur = 1
                    else:
                        keepInd.append(j)
                mergedCluster = [ mergedCluster[i] for i in keepInd ]
#                print "merged cluster length: " + str(len(mergedCluster))
        sizeMergedCluster = np.zeros(len(mergedCluster))
        for i in range(0,len(mergedCluster)):
            sizeMergedCluster[i] = len(mergedCluster[i])
        sortSize = np.sort(sizeMergedCluster)[::-1]
        sortMergedInd = np.argsort(sizeMergedCluster)[::-1]
        mergedCluster = [ mergedCluster[i] for i in sortMergedInd ]
        currentInd = 0
    return mergedCluster
        

def mainroutine(step1, inputdata, nRow, nCol, gamma, t, lambda_val, beta, minClusterSize, method = "Pearson"):
    if step1:
        # step 1: get PCC mat
        tm = time.time()
        cMatrix = massivePCC(inputdata, nRow, nCol, method)
        print "step 1: get " + method + " Correlation Cefficient mat elapsed time: " + str("{0:.2f}".format(time.time() - tm)) + " seconds"
    else:
        cMatrix = inputdata
    
    # step 2: calculate lmQCM
    tm = time.time()
    cell = localMaximumQCM(abs(cMatrix), gamma, t, lambda_val)
    print "step 2: calculate lmQCM elapsed time: " + str("{0:.2f}".format(time.time() - tm)) + " seconds"
    
    # step 3: Merge the overlapped networks
    tm = time.time()
    mergedCluster = lmQCMMerge(cell, beta, minClusterSize)
    print "step 3: Merge the overlapped networks elapsed time: " + str("{0:.2f}".format(time.time() - tm)) + " seconds"
    print "Length of merged cluster: " + str(len(mergedCluster))
    return mergedCluster


#import pickle
#with open('/Users/zhi/Desktop/GeneCoexpression/shiny/finalExpArray.pkl', 'rb') as f:
#    finalExpArray = pickle.load(f)
#gamma = 0.5
#t = 1
#lambda_val = 1
#beta = 0.4
#minClusterSize = 10
#mergedCluster = mainroutine(1,finalExpArray, 18172, 415, gamma, t, lambda_val, beta, minClusterSize)