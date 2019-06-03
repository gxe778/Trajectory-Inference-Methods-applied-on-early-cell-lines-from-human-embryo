###############################################################################
##
## Demonstration script
##
## SLICE Analysis of HEE data
##
## SLICE version: SLICE a12242016
##
## Author: Minzhe Guo (minzhe.guo@cchmc.org)
## Date: Dec. 24, 2016
##
###############################################################################



# loading SLICE functions
source("slice.R")


# use a variable to store the path to the data directory; getwd() functions returns the full path of current working directory
data.dir <- paste(getwd(),"/data/", sep="")

# use a variable to store the name of the dataset and analysis
data.name <- "HEE"


cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Start SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))


# a string describing the analysis and will be attached to the names of all the files generated during the execution
context_str = paste("SLICE-", data.name, "-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="") 


###########################################################################
#####                                                               #######
#####               Loading data and preprocessing                  #######
#####                                                               #######
###########################################################################


# loading the HEE data set; expression profiles with same gene symbol annotations have been averaged
load(paste(data.dir, "HEE.Rda", sep=""))

# the expressions, cell info, and gene info are stored as an Biobase::ExpressionSet object
es <- HEE

# remove two morula cells that were consistently, and can be considered as outliers based on PCA (Yan et al., 2013) and hierarchical clustering (Guo et al., 2015)
es <- es[, which(!(rownames(pData(es)) %in% c("Morulae #1 -Cell#3", "Morulae #1 -Cell#8")))]

# update the factor levels of pData(es)$Cluster: human early embryo developmental stage information of cells
pData(es)$Cluster <- factor(pData(es)$Cluster)


# remove ERCC genes and Ribosomal genes
ercc.genes <- grep("^ERCC-", rownames(fData(es)), value = TRUE)
rb.genes <- grep("^RPL|^RPS|^MRPL|^MRPS", rownames(fData(es)), value = TRUE)
es <- es[which(!(rownames(fData(es)) %in% c(rb.genes, ercc.genes))), ]


###########################################################################
#####                                                               #######
#####                   Initlializing SLICE object                  #######
#####                                                               #######
###########################################################################

# use the expression matrix of HEE cells (n=88) to create a SLICE object
# use the development stage information from the original analysis (encoded in the Cluster attribute of the pData) as the original identity
# set projname to describe this analysis
sc <- construct(exprmatrix=as.data.frame(exprs(es)), 
                cellidentity=factor(pData(es)$Cluster, levels=c("Oocyte","Zygote","2-cell embry","4-cell embry", "8-cell embry","Morulae","Late blastocyst")),
                projname=context_str)


###########################################################################
#####                                                               #######
#####     Measuring cell differentiation state using scEntropy      #######
#####                                                               #######
###########################################################################


# loading the pre-computed human gene-gene Kappa similarity matrix
# The matrix will be loaded into a variable named "km"
load(paste(data.dir, "hs_km.Rda", sep=""))


# bootstrap estimation of scEntropy. 100 boostrapped samples, 1000 genes in each sample
sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of human genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=100,                             # 100 iterations
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=201602)                    # set the random seed to reproduce the results in the paper

# plot the entropy; will be save to a pdf file named with "entropies" suffix in the working directory
plotEntropies(sc)


###########################################################################
#####                                                               #######
#####                      Lineage reconstruction                   #######
#####                                                               #######
###########################################################################


# perform PCA using the expression of genes expressed in at least 30% of cells and with greater than 0.5 variance in log2 RPKM
sc <- getRDS(sc, method="pca", num_dim=2, log.base=2, do.center=TRUE, do.scale=FALSE, use.cor=F, min.var=0.5, min.cells=NULL)


# use the graph-based method in SLICE to infer lineage model
# results will be plotted on screen and saved in the PDF files with file names containing "lineage-g"
sc <- getLineageModel(sc, lm.method="graph",                          # select graph-based method
                      model.type="tree",                              # infer mst based lineage model
                      reverse=F,                                      # infer differentiation lineage
                      ss.method="pcst", ss.threshold=0.25,            # use linear prize collecting tree algorithm to find core cell sets
                      wiring.threshold=function(mst) max(mst)*0.2,    # the threshold for local wiring
                      community.method="louvain")                     # use the louvain algorithm to detect cell clusters

# marker genes for lineage validate
markers <- c("DPPA3", "SALL2", "KLF5", "LIN28A", "LIN28B", "FGFR4", "CLDN10")

# use the shortest-path method in SLICE to reconstruct cell transitional path following C1->C2->C3->C4
# results will be plotted on screen and saved in the PDF files with file names containing "path-sp"
sc <- getTrajectories(sc, method="sp", start=1, end=4, network="mst", NN.threshold=0.8, NN.type="all", do.plot=T)

# extract and visualize expression profiles of marker genes in the reconstructed shrotest-path transitional path
# results will be plotted on screen and saved in the PDF files with file names containing "path-sp" and "profiles"
sc <- getProfiles(sc, trajectory.type="sp", genes.use=markers)


# use the clustering-based method in SLICE to infer lineage model
# results will be plotted on screen and saved in the PDF files with file names containing "lineage-c"
sc <- getLineageModel(sc, lm.method="clustering", model.type="tree", reverse=F, ss.threshold=0.25, cluster.method="pam", k=4)

# use the principal-curve method in SLICE to reconstruct cell transitional path following C1->C2->C3->C4
# results will be plotted on screen and saved in the PDF files with file names containing "path-pc"
sc <- getTrajectories(sc, method="pc", start=1, end=4, do.trim=F,  do.plot=T)

# extract and visualize expression profiles of marker genes in the reconstructed principal-curve transitional path
# results will be plotted on screen and saved in the PDF files with file names containing "path-pc" and "profiles"
sc <- getProfiles(sc, trajectory.type="pc", genes.use=markers)


cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Finish SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))