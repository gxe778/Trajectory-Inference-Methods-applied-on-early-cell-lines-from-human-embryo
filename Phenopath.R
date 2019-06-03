#load HEE
load("HEE.Rda")

#load PhenoPath
install.packages("devtools")
devtools::install_github("kieranrcampbell/phenopath", build_vignettes = TRUE)
library(phenopath)

source("https://bioconductor.org/biocLite.R")
biocLite(c("scater", "MultiAssayExperiment"))
library(MultiAssayExperiment)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)



#create se
se <- SingleCellExperiment(assays=list(count=HEE@assayData$exprs)) 

as(se,"SingleCellExperiment")
library(phenopath)

x<-colnames(se)
x<-factor(x, levels=c(unique(x)))
scale_vec <- function(x) (x - mean(x)) / sd(x)
pc1 <- prcomp(t(HEE@assayData$exprs))$x[,1]
pc1 <- scale_vec(pc1)
time_numeric <- (gsub("h", "", x))
pc1 <- pc1 * sign(cor(pc1, x))

#determine trajectory
y <- as.matrix(HEE@assayData$exprs)
fit <- phenopath(y, x, z_init = pc1)
trajectory(fit)
