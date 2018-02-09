# 01/31/2018 Zhi Huang
# GEOmetadb
library(GEOmetadb)

if(!file.exists('GEOmetadb.sqlite')) {
  getSQLiteFile()
}

file.info('GEOmetadb.sqlite')
