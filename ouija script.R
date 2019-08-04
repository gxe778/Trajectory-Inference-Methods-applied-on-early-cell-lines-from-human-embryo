#Ouija Script

devtools::install_github("kieranrcampbell/ouija")
library(ouija)



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


data(v) # synthetic gene expression data bundled
oui <- ouija(v)
pseudotimes <- map_pseudotime(v)
