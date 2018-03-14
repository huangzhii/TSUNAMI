import pickle
def mainroutine(step1, inputdata, nRow, nCol, gamma, t, lambda_val, beta, minClusterSize, method = "Pearson"):
  pickle.dump( inputdata, open( "/Users/zhi/Desktop/inputdata.pkl", "wb" ) )
