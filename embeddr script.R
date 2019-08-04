#Embeddr R script
#https://github.com/kieranrcampbell/embeddr

install_github('kieranrcampbell/embeddr')

library(embeddr)

#https://github.com/kieranrcampbell/embeddr

library(embeddr)

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

# create a SingleCellExperiment object sce
sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))



## Create a cell-cell correlation graph and use it for the reduced embedding:
sce <- embeddr(sce)

sce<-laplacian_eigenmap(v)

## Plot a reduced-dimension embedding
plot_embedding(v)

## Optionally cluster the embedding. Cluster assignments are stored in pData(sce)$cluster.
## If no number of clusters is designated, the number is chosen using the BIC from package mclust
sce <- cluster_embedding(sce)

## Fit pseudotime using principal curves
sce <- fit_pseudotime(sce)

## Plot genes 1:10 in pseudotime:
plot_in_pseudotime(sce[1:10,])

## Fit differential expression pseudotime model. This will report gene name, p-val and q-val
diff_gene_test <- pseudotime_test(sce)
