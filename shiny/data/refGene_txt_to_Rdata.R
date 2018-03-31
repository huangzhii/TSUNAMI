# 03/30/2018 Zhi Huang
# convert refGene data from raw txt to Rdata
# original data downloaded from:
# HG19: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
# HG38: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

hg19 <- read.csv("./data/hg19refGene_20180330.txt",sep = "\t", header=F)
hg38 <- read.csv("./data/hg38refGene_20180330.txt",sep = "\t", header=F)

save(hg19,file="UCSC_hg19_refGene_20180330.Rdata")
save(hg38,file="UCSC_hg38_refGene_20180330.Rdata")
