                                              #phenopath 


#devtools::install_github("kieranrcampbell/phenopath", build_vignettes = TRUE)

library(phenopath)



#load data

TI_Data_Yan_et_al <- read_excel("TI_Yan.xlsx")



TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)





TID<-TI_Data_Yan_et_al





rownames(TID)<-TI_Data_Yan_et_al[,1]





TID<-TID[,-c(1:2)]





#save data in excel and R:


    save(TID,file="TID.Rda")



    write.xlsx2(TID, file="TID.xlsx")




counts<-as.matrix(TID)



v<-log2(counts+1)


#save: excel and R:


   save(v,file="v.Rda")


    write.xlsx2(v, file="v.xlsx")



#create a SingleCellExperiment object sce


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))


#pre-processing

x <- 2 * (v== "v") - 1



scale_vec <- function(x) (x - mean(x)) / sd(x)

pc1 <- prcomp(t(x))





#phenopath

fit <- phenopath(v, x)
