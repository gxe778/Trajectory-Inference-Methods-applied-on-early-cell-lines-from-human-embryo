#CellTrails: HEE 

#install CellTrails from bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("CellTrails", version = "3.8")

#load CellTrails
library(CellTrails)

#cellTrails uses Bioconductor single cell experiment class (Lun and Risso 2017): install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)

#load HEE data:
setwd("C:/Users/G177/Desktop/Trajectory Inference Repetoire")
load("HEE.Rda")


#create single cell experiment with HEE data(in S4 Biobase format)
HEESCE<-SingleCellExperiment(assays=list(logcounts=HEE@assayData$exprs), colData=HEE@phenoData@data)

#Use this function repeatedly to check tracctory inference of data:
showTrajInfo(HEESCE)

#filter out top 100 genes with maximal variance
library(magrittr)
assay(HEESCE) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(HEESCE)
vars <- sort(vars, decreasing = TRUE)

HEESCE <- HEESCE[names(vars[1:9000]),]

#filter by detection level
tfeat <- filterTrajFeaturesByDL(HEESCE, threshold=0)

#set features to objects
trajFeatureNames(HEESCE) <- tfeat
showTrajInfo(HEESCE)


#use spectral embedding to create weighted graph of cells (i.e. manifold: nodes represent cells and edges connecting nodes represent high statisitcal dependency(i.e. similiar in expression profile)):
se <- embedSamples(HEESCE)
names(se)

#determine dimensions of manifold in higher space using all 90 of eigenvalues in se:
d <- findSpectrum(se$eigenvalues, frac=90)

# set latent space to object:
latentSpace(HEESCE) <- se$components[, d]


# obtain manifold
gp<-plotManifold(HEESCE, color_by="phenoName", name="Sample", recalculate = TRUE)
manifold2D(HEESCE) <- gp
gp


#determine clusters
cl <- findStates(HEESCE, min_size=.01, min_feat=5, max_pval=1e-4, min_fc=2)
head(cl)
states(HEESCE)<-cl

#plot manifold with clusters:
plotManifold(HEESCE, color_by="phenoName", name="state")
plotStateSize(HEESCE)

#traj. graph 
HEESCE<- connectStates(HEESCE, l=20)

#component determination(i.e. certain pathway)
trajComponents(HEESCE)


gp<-plotStateTrajectory(HEESCE, color_by="phenoName", name="Cluster",
                    component=1,recalculate=TRUE)
gp

#trajectory inference for certain component(3 total)
HEESCE <- selectTrajectory(HEESCE, component=1)

HEESCE_subset <- HEESCE[, trajSampleNames(HEESCE)]

plotStateSize(HEESCE_subset)

HEESCE <- fitTrajectory(HEESCE)
showTrajInfo(HEESCE)

plotTrajectoryFit(HEESCE)

#graph map of inference
write.ygraphml(sce=HEESCE,
               file='HEESCE.graphml',
               color_by='phenoName',
               name='state',
               node_label='state')

#use yED live to construct graph

#save map
tl <- read.ygraphml("HEESCE.graphml")

#plot layout
plot(tl[,1:2], axes=FALSE, xlab="", ylab="", pch=20, cex=.25)

trajLayout(HEESCE, adjust=TRUE) <- tl

showTrajInfo(HEESCE)

plot(trajLayout(HEESCE), axes=FALSE, xlab="", ylab="", pch=20, cex=.25)

plotMap(HEESCE, color_by="phenoName", name="Sample")

#use gene for topology: CCND2
plotMap(HEESCE, color_by="featureName", name="CCND2", type="surface.fit")

#trail definition landmarks
HEESCE <- addTrail(HEESCE, from="B1", to="H1", name="Tr1")
HEESCE <- addTrail(HEESCE, from="B1", to="H2", name="Tr2")

plotTrail(HEESCE, name="Tr1")
plotTrail(HEESCE, name="Tr2")

plotManifold(HEESCE, color_by="phenoName", name="Tr1")

#pseudotime
ptime <- trails(HEESCE)[, "Tr1"]

#subset samples
Tr1 <- HEESCE[, !is.na(ptime)]

o <- order(trails(tr1)[, "Tr1"])
Tr1 <- Tr1[, o]
ptime <- trails(Tr1)[, "Tr1"]
names(ptime) <- colnames(Tr1)

#ptime per state
ptime_states <- split(ptime, states(Tr1))
lptime <- lapply(ptime_states, 
                 function(x){y <- diff(sort(x)); y[-length(y)]})

#trail identification
plotMap(HEESCE, color_by="phenoName", name="landmark")

plotTrail(HEESCE, name="Tr1")
