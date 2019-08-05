#install functions 


install.packages("reticulate")
library(reticulate)

install.packages("dplyr")
library(dplyr)

install.packages("devtools")
library(devtools)

install.packages("githubinstall")
library(githubinstall)

devtools::install_github("xu-lab/SINCERA")
library(SINCERA)

install.packages("foreign")
library(foreign)

install.packages("ggplot2")
library(ggplot2)

install.packages("MASS")
library(MASS)

install.packages("boot")
library(boot)

install.packages("cluster")
library(cluster)

install.packages("rgl")
library(rgl)

install.packages("VGAM")
library(VGAM)

install.packages("plotCluster")
library(plotCluster)

install.packages("DDRTree")
library(DDRTree)

install.packages("covr")
library(covr)

install.packages("gam")
library(gam)

install.packages("KNNShiny")
library(KNNShiny)

install.packages("pheatmap")
library(pheatmap)

install.packages("fastICA")
library(fastICA)

install.packages("densityClust")
library(densityClust)

install.packages("mclust")
library(mclust)

install.packages("reshape")
library(reshape)

install.packages("lubridate")
library(lubridate)

install.packages("Rdimtools")
library(Rdimtools)

install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("tidyverse")
library(tidyverse)

install.packages("tidyr")
library(tidyr)

install.packages("shiny")
library(shiny)

install.packages("plyr")
library(plyr)

install.packages("graph")
library(graph)

install.packages("Seurat")
library(Seurat)


#Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install()



BiocManager::install("graph")
library(graph)

BiocManager::install("limma")
library(limma)

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("annotation")
library(annotation)

BiocManager::install("sequencing")
library(sequencing)

BiocManager::install("arrays")
library(arrays)

BiocManager::install("ExpressionNormalizationWorkflow")
library(ExpressionNormalizationWorkflow)

BiocManager::install("simpleSingleCell")
library(simpleSingleCell)

BiocManager::install("cytofWorkflow")
library(cytoworkflow)

BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

BiocManager::install("SCnorm")
library(SCnorm)

BiocManager::install("scruff")
library(scruff)

BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

BiocManager::install("Rgraphviz")
library(Rgraphviz)

BiocManager::install("ensembldb")
library(ensembldb)

BiocManager::install("Geoquery")
library(Geoquery)

BiocManager::install("pcaMethods")
library(pcaMethods)

BiocManager::install("KEGGgraph")
library(KEGGgraph)

BiocManager::install("scater")
library(scater)

BiocManager::install("destiny")
library(destiny)

BiocManager::install("metagenomeSeq")
library(metagenomeSeq)

BiocManager::install("flowWorkspace")
library(flowWorkspace)

BiocManager::install("zinbwave")
library(zinbwave)

BiocManager::install("arrayQualityMetrics")
library(arrayQualityMetrics)

BiocManager::install("flowclust")
library(flowclust)

BiocManager::install("clusterExperiment")
library(clusterExperiment)

BiocManager::install("genomation")
library(genomation)

BiocManager::install("GenVisR")
library(GenVisR)

BiocManager::install("Bionet")
library(BioNet)


BiocManager::install("SeqGSEA")
library(SeqGSEA)

BiocManager::install("KEGGprofile")
library(KEGGprofile)

BiocManager::install("scmap")
library(scmap)


BiocManager::install("PharmacoGx")
library(PharmacoGx)

BiocManager::install("BASiCS")
library(BASiCS)

BiocManager::install("GenomicInteractions")
library(GenomicInteractions)

BiocManager::install("GOSim")
library(GOSim)

BiocManager::install("AnnotationHubData")
library(AnnotationHubData)


BiocManager::install("BiocWorkflowTools")
library(BiocWorkflowTools)


BiocManager::install("ClusterSignificance")
library(ClusterSignificance)

BiocManager::install("diffGeneAnalysis")
library(diffGeneAnalysis)

BiocManager::install("GSCA")
library(GSCA)

BiocManager::install("CellMapper")
library(CellMapper)

BiocManager::install("GSReg")
library(GSReg)

BiocManager::install("SeqSQC")
library(SeqSQC)

BiocManager::install("scFeatureFilter")
library(scFeatureFilter)

BiocManager::install("GIGSEA")
library(GIGSEA)

BiocManager::install("PCAtools")
library(PCAtools)

BiocManager::install("batchelor")
library(batchelor)

BiocManager::install("pathwayPCA")
library(pathwayPCA)

BiocManager::install("CellBench")
library(CellBench)

BiocManager::install("bayNorm")
library(bayNorm)

BiocManager::install("scFeatureFilter")
library(scFeatureFilter)

BiocManager::install("HCABrowser")
library(HCABrowser)

BiocManager::install("scRecover")
library(scRecover)

BiocManager::install("scAlign")
library(scAlign)

BiocManager::install("bigPint")
library(bigPint)

BiocManager::install("splatter")
library(splatter)

BiocManager::install("BioNet")
library(BioNet)

BiocManager::install("sincell")
library(sincell)

BiocManager::install("scone")
library(scone)
