# from __future__ import division
import pickle
import numpy as np
import json

def massivePCC(finalExpArray, nRow, nCol):

	# with open('./PCCconsise.pkl', 'rb') as f:
	# 	PCCconsise = pickle.load(f)
	# return PCCconsise.tolist()
	finalExp = np.asarray(finalExpArray).reshape((nCol, nRow)).T
	if nRow > 1:
		sumX = np.sum(finalExp, axis=1)
		sumV = np.sqrt( np.sum(np.multiply(finalExp,finalExp), axis=1) - np.multiply(sumX,sumX)/nCol )
		PCC_mat = np.zeros((nRow, nRow))
		for i in range(0,nRow-1):
			if (i % (nRow//10)) == 0:
				percentile = 100*i/nRow
				print "Processing Massive PCC ... " + str(percentile) + "%"
			temp = np.divide(np.multiply(sumX[i], sumX[(i+1):nRow]), nCol)
			np.divide(( np.dot(finalExp[(i+1):nRow,:], finalExp[i,:]) - temp ) \
				, np.multiply(sumV[i],sumV[(i+1):nRow]), out = PCC_mat[i,(i+1):nRow])
			# PCC_mat[i,(i+1):nRow] = ( np.dot(finalExp[(i+1):nRow,:], finalExp[i,:]) - np.multiply(sumX[i], sumX[(i+1):nRow])/nCol ) / np.multiply(sumV[i],sumV[(i+1):nRow])
		np.add(PCC_mat, PCC_mat.T, out=PCC_mat)

	return PCC_mat.tolist()

# with open('./finalExpArray.pkl', 'rb') as f:
# 	finalExpArray = pickle.load(f)
# finalExpArray = 0
# massivePCC(finalExpArray, 18172, 415)