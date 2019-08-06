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




#Load Data
TI_Data_Yan_et_al <- read_excel("TI Data Yan et al.xlsx")


TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)


TID<-TI_Data_Yan_et_al


rownames(TID)<-TI_Data_Yan_et_al[,1]


TID<-TID[,-c(1:2)]


#save data in excel and R
save(TID, file="TID.xlsx")

save(TID,file="TID.Rda")


#Create SingleCell Object
counts<-as.matrix(TID)

v<-log2(counts+1)

sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))

#slice functions
#make directory SLICE

source("slice.R")




# use a variable to store the name of the dataset and analysis
data.name <- "V.Rda"


cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Start SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))


# a string describing the analysis and will be attached to the names of all the files generated during the execution
context_str = paste("SLICE-", data.name, "-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="") 





###########################################################################
#####                                                               #######
#####                   Initlializing SLICE object                  #######
#####                                                               #######
###########################################################################

# use the expression matrix of HEE cells (n=88) to create a SLICE object
# use the development stage information from the original analysis (encoded in the Cluster attribute of the pData) as the original identity
# set projname to describe this analysis
sc <- construct(exprmatrix=as.data.frame(v), 
                cellidentity=factor(colnames(v)),
                projname=context_str)


###########################################################################
#####                                                               #######
#####     Measuring cell differentiation state using scEntropy      #######
#####                                                               #######
###########################################################################


# loading the pre-computed human gene-gene Kappa similarity matrix
# The matrix will be loaded into a variable named "km"
 load("hs_km.Rda")


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



# use the shortest-path method in SLICE to reconstruct cell transitional path: ex. C2:C4
sc <- getTrajectories(sc, method="sp", start=2, end=4, network="mst", NN.threshold=0.8, NN.type="all", do.plot=T)



# use the clustering-based method in SLICE to infer lineage model
sc <- getLineageModel(sc, lm.method="clustering", model.type="tree", reverse=F, ss.threshold=0.25, cluster.method="pam", k=6)

# use the principal-curve method in SLICE to reconstruct cell transitional path: ex. C5:C1
sc <- getTrajectories(sc, method='sp', start=5, end=1, do.trim=F,  do.plot=T)




cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Finish SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))
