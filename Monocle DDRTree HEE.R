#HEE: Monocle

load ("HEE.Rda")

source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("monocle")

devtools::install_github("cole-trapnell-lab/DDRTree", ref="simple-ppt-like")

devtools::install_github("cole-trapnell-lab/DDRTree", ref="simple-ppt-like")

devtools::install_github("cole-trapnell-lab/L1-graph")


install.packages("reticulate")
library(reticulate)
py_install('umap-learn', pip = T, pip_ignore_installed = T) 
# Ensure the latest version of UMAP is installed
py_install("louvain")

devtools::install_github("cole-trapnell-lab/monocle-release", ref="monocle3_alpha")

library(monocle)

row.has.na <- apply(HEE, 1, function(x){any(is.na(x))})

sum(row.has.na)

HEE <- HEE[!row.has.na,]


HEE<-newCellDataSet(HEE@assayData$exprs, phenoData = HEE@phenoData, featureData = HEE@featureData)

HEE <- estimateSizeFactors(HEE)
HEE <- estimateDispersions(HEE)


HEE <- preprocessCDS(HEE, num_dim = 20)


HEE <- DDRTree(HEE@assayData$exprs)

HEE <- partitionCells(HEE)

cds <- learnGraph(HEE,  RGE_method = "DDRTree")


celltype<-HEE@phenoData@data$Sample

plot_cell_trajectory(HEE,
                     color_by = "cell_type") +
  scale_color_manual(values = cell_type_color)
