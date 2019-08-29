#Embeddr R script

#https://github.com/kieranrcampbell/embeddr



install_github('kieranrcampbell/embeddr')



library(embeddr)



#Load Data


TI_Data_Yan_et_al <- read_excel("TI Data Yan et al.xlsx")




TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)





TID<-TI_Data_Yan_et_al





rownames(TID)<-TI_Data_Yan_et_al[,1]




TID<-TID[,-c(1:2)]




#save data in excel and R

save(TID,file="TID.Rda")



write.xlsx2(TID, file="TID.xlsx")





#create a SingleCellExperiment object sce



counts<-as.matrix(TID)



v<-log2(counts+1)





sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))




#Create a cell-cell correlation graph and use it for the reduced embedding:

sce <- embeddr(sce)



## Plot a reduced-dimension embedding

plot_embedding(sce)



#clustering

sce <- cluster_embedding(sce)



## Fit pseudotime using principal curves

sce <- fit_pseudotime(sce)





## Fit differential expression pseudotime model. This will report gene name, p-val and q-val

diff_gene_test <- pseudotime_test(sce)

