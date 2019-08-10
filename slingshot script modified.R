#Slingshot script

#https://github.com/kstreet13/slingshot

BiocManager::install("kstreet13/slingshot")

library(slingshot)

#Load Data
TI_Data_Yan_et_al <- read_excel("TI Data Yan et al.xlsx")


TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)


TID<-TI_Data_Yan_et_al


rownames(TID)<-TI_Data_Yan_et_al[,1]


TID<-TID[,-c(1:2)]


#save data in excel and R
save(TID,file="TID.Rda")

write.xlsx2(TID, file="TID.xlsx")


# create a SingleCellExperiment object sce

counts<-as.matrix(TID)

v<-log2(counts+1)


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))



#Dimension Reduction

PCA <- prcomp(t(log1p(v)))


#Slingshot
sim <- slingshot(v,  reducedDim='PCA')
