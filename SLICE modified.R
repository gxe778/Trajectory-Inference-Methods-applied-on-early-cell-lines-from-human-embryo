                              #SLICE



#load data

TI_Data_Yan_et_al <- read_excel("TI_Yan.xlsx")




TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)





TID<-TI_Data_Yan_et_al





rownames(TID)<-TI_Data_Yan_et_al[,1]





TID<-TID[,-c(1:2)]





#save data in excel and R:


            #save(TID,file="TID.Rda")



            #write.xlsx2(TID, file="TID.xlsx")






counts<-as.matrix(TID)



v<-log2(counts+1)


#save data in excel and R:


            #save(v, file="v.Rda")



            #write.xlsx2(v, file="v.xlsx")




#create SingleCell object



sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))





#download SLICE: #https://research.cchmc.org/pbge/slice.html





#SLICE functions: "slice.R": open and source




#use a variable to store the name of the dataset and analysis

data.name <- "v.Rda"





cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Start SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))





#a string describing the analysis and will be attached to the names of all the files generated during the execution

context_str = paste("SLICE-", data.name, "-", format(Sys.time(), "vSLICE"), "-", sep="") 




###########################################################################

#####                                                               #######

#####                   Initlializing SLICE object                  #######

#####                                                               #######

###########################################################################



#use the expression matrix (v) to create a SLICE object

sc <- construct(exprmatrix=as.data.frame(v), 
                
                cellidentity=factor(colnames(v)),
                
                projname=context_str)





###########################################################################

#####                                                               #######

#####     Measuring cell differentiation state using scEntropy      #######

#####                                                               #######

###########################################################################





#load pre-computed human gene-gene Kappa similarity matrix: "km"

load("hs_km.Rda")



#bootstrap estimation of scEntropy: 100 boostrapped samples, 1000 genes in each sample

sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of human genes
                 
                 calculation="bootstrap",               # choose the bootstrap calculation
                 
                 B.num=100,                             # 100 iterations
                 
                 exp.cutoff=1,                          # the threshold for expressed genes
                 
                 B.size=1000,                           # the size of bootstrap sample
                 
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 
                 random.seed=201602)                   # set the random seed to reproduce the results in the paper



#plot entropy; will save to a pdf file named with "entropies" suffix in the working directory

plotEntropies(sc)





###########################################################################

#####                                                               #######

#####                      Lineage reconstruction                   #######

#####                                                               #######

###########################################################################





#perform PCA using the expression of genes expressed in at least 30% of cells and with greater than 0.5 variance in log2 RPKM

sc <- getRDS(sc, method="pca", num_dim=2, log.base=2, do.center=TRUE, do.scale=FALSE, use.cor=F, min.var=0.5, min.cells=NULL)





#graph-based method in SLICE to infer lineage model; results plotted on screen and saved in PDF files: "lineage-g"


sc <- getLineageModel(sc, lm.method="graph",                          #select graph-based method
                      
                      model.type="tree",                              #infer mst based lineage model
                      
                      reverse=F,                                      #infer differentiation lineage
                      
                      ss.method="pcst", ss.threshold=0.25,            #use linear prize collecting tree algorithm to find core cell sets
                      
                      wiring.threshold=function(mst) max(mst)*0.2,    #threshold for local wiring
                      
                      community.method="louvain")                     #louvain algorithm to detect cell clusters







#shortest-path method in SLICE to reconstruct cell transitional path

sc <- getTrajectories(sc, method="sp", start=2, end=4, network="mst", NN.threshold=0.8, NN.type="all", do.plot=T)







#clustering-based method in SLICE to infer lineage model

sc <- getLineageModel(sc, lm.method="clustering", model.type="tree", reverse=F, ss.threshold=0.25, cluster.method="pam", k=6)



#principal-curve method in SLICE to reconstruct cell transitional path

sc <- getTrajectories(sc, method='sp', start=5, end=1, do.trim=F,  do.plot=T)









cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Finish SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))



