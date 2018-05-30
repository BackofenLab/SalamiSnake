source('http://bioconductor.org/biocLite.R')
biocLite('RCAS')

library(RCAS)

args = commandArgs(trailingOnly = TRUE)
runReport(
queryFilePath = args[1],
gffFilePath = args[2],
#msigdbFilePath = args[3], if msigdbAnalysis==TRUE
annotationSummary = TRUE,
goAnalysis = TRUE,
msigdbAnalysis = FALSE,
motifAnalysis = FALSE,
genomeVersion = args[4],
outDir = args[5],
printProcessedTables = TRUE,
sampleN = args[6],
selfContained = TRUE)
