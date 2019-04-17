# source('http://bioconductor.org/biocLite.R')
# biocLite('RCAS')

source("https://bioconductor.org/biocLite.R")
biocLite("org.Dm.eg.db")

args = commandArgs(trailingOnly=TRUE)

library(RCAS)
library(org.Dm.eg.db)

print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

args = commandArgs(trailingOnly = TRUE)
runReport(
queryFilePath = args[1],
gffFilePath = args[2],
#msigdbFilePath = args[3], if msigdbAnalysis==TRUE
annotationSummary = TRUE,
goAnalysis = TRUE,
msigdbAnalysis = FALSE,
motifAnalysis = FALSE,
genomeVersion = args[3],
outDir = args[4],
printProcessedTables = TRUE,
sampleN = args[5],
selfContained = TRUE)
