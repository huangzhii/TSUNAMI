# A sample .Rprofile file with two different package repositories.
local({
  r <- getOption("repos")
  r["bioconductor"] <- "https://www.bioconductor.org/"
  r["CRAN"] <- "https://cran.rstudio.com/"
  r["CRAN2"] <- "http://ftp.ussg.iu.edu/CRAN/"
  options(repos = r)
  # options(useHTTPS=FALSE, BioC_mirror="http://bioconductor.org")
  # source("http://bioconductor.org/biocLite.R")
})

