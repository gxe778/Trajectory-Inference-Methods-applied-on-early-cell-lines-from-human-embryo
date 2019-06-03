#EmbeddR: HEE dataset

#load HEE dataset:
load("HEE.Rda")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


install_github('davismcc/scater')
BiocManager::install('devtools')
BiocManager::install("SingleCellExperiment")
install_github('kieranrcampbell/embeddr')

library(scater)
library(SingleCellExperiment)
library(embeddr)
library(devtools)




#create sce
pd <- new('AnnotatedDataFrame', data=HEE@phenoData@data)

sce <- SingleCellExperiment(assays=HEE@assayData$exprs, colData=HEE@phenoData@data)


sce <- fit_pseudotime(sce)

diff_gene_test <- pseudotime_test(sce)

library(embeddr)
sce <- embeddr(sce)
