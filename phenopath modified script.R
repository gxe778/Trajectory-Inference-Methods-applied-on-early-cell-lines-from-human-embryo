#phenopath script

#https://github.com/kieranrcampbell/phenopath

devtools::install_github("kieranrcampbell/phenopath", build_vignettes = TRUE)

library(phenopath)

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


#Phenopath pre-processing
x <- 2 * (v== "v") - 1

scale_vec <- function(x) (x - mean(x)) / sd(x)
pc1 <- prcomp(t(x))


#Phenopath
fit <- phenopath(v, x)
