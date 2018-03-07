# test1.R
Sys.setenv(PATH = paste("/Users/zhi/anaconda2/bin", Sys.getenv("PATH"), sep=":")) 
Sys.setenv(PATH = paste("usr/local/bin", Sys.getenv("PATH"), sep=":")) 

library(rPython)

setwd("~/Desktop")
system("which python")
system("python --version")
python.load("test1.py")
system("/Users/zhi/anaconda2/bin/python test1.py")
system("/usr/local/bin/python test1.py")
