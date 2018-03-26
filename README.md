# TSUNAMI: Translational Bioinformatics Tool SUite for Network Analysis and MIning
Gene co-expression network (GCN) mining aims to mine gene modules with highly correlated expression profiles across sample cohorts. It may help to reveal latent molecular interactions, identify novel gene functions, pathways and drug targets, as well as providing disease mechanistic insights on for biological researchers. TSUNAMI is developed to allow biological researchers with no programing background to perform GCN mining themselves. Users can get access to our online service from https://apps.medgen.iupui.edu/rsc/tsunami/

## Features
It has several highlight features and advantages:
* User friendly interface, easy-access and real-time co-expression network mining based on web server;
* Direct access and search of GEO database as well as user-input expression matrix for network mining;
* Support multiple data formats and data preprocessing interface is bundled together;
* Multiple co-expression analysis tools available with a high flexibility of variable selection;
* Integrated downstream Enrichr GO enrichment analysis and link to other GO tools as well;
* All results can be downloaded with multiple formats (CSV, txt, etc.).

## Required R packages
* shiny
* rsconnect
* plyr
* data.table
* genefilter
* Biobase
* lmQCM
* WGCNA
* GEOquery
* dplyr
* enrichR
* DT
* plotly
* openxlsx
* survival

## Installation
* Install R and RStudio

Please download the right version of R for your own system from https://cloud.r-project.org/ and install it. 
RStudio is an integrated development environment (IDE) for R. Please download RStudio from https://www.rstudio.com/products/rstudio/#Desktop and install it.
* Install R packages within RStudio as below.
```r
install.packages(c("shiny", "rsconnect", "plyr", "data.table", "lmQCM", "WGCNA", "dplyr", "enrichR", "DT", "plotly", "openxlsx", "survival"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("genefilter", "Biobase", "AnnotationDbi", "impute", "GO.db", "preprocessCore", "GEOquery"))
```
* Download the source code for TSUNAMI from GitHub https://github.com/huangzhii/TSUNAMI/.

Open ui.R from the download TSUNAMI folder. Click the button named "Run App" on the upper right corner of code editing window, an web page will be automatically invoked and the TSUNAMI application is ready to use. It may take a few minutes to load the required packages. 
