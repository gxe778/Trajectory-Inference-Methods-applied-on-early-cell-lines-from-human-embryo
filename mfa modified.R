#mfa Script
#https://github.com/kieranrcampbell/mfa

#Install
devtools::install_github("kieranrcampbell/mfa", build_vignettes = TRUE)

library(mfa)

#Load Data Yan et al.
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

#Run mfa
m <- mfa(v)
