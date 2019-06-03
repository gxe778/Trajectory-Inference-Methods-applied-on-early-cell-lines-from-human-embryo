#DPT: Diffusion Pseudotime: HEE: Thtrough R package Destiny

#Download Diffusion Mapping
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")library


BiocManager::install("destiny")
library (destiny)

#load HEE
load("HEE.Rda")
data("guo")

#Create DPT Object
dm <- DiffusionMap(HEE) 
dpt <- DPT(dm)
plot(dpt)

#cluster and draw path from 1 to 4:
plot(dpt, root = 2, paths_to = c(1,4), col_by = 'branch')
