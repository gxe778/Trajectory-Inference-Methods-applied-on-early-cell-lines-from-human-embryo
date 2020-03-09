                # tSpace

                # https://github.com/hylasD/tSpace





# Install/Load tSpace


#devtools::install_github('hylasD/tSpace', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), force = T)


library(tSpace)

#Load Data

TI_Data_Yan_et_al <- read_excel("TI_Yan.xlsx")



TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)





TID<-TI_Data_Yan_et_al





rownames(TID)<-TI_Data_Yan_et_al[,1]





TID<-TID[,-c(1:2)]





# save data in excel and R:


# save(TID,file="TID.Rda")



# write.xlsx2(TID, file="TID.xlsx")




counts<-as.matrix(TID)



v<-log2(counts+1)


# save: excel and R:


# save(v,file="v.Rda")


# write.xlsx2(v, file="v.xlsx")



# create a SingleCellExperiment object sce


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))



# run tSpace 


ts <- tSpace(df = v, K = 20, L = 15, D = 'pearson_correlation', graph = 5, trajectories = 200, wp = 15, dr = 'pca', core_no = 2)


# tSpace: 1000 trajectories (i.e. cell # < 10000) 


ts <- tSpace(df = v, K = 20, L = 15, D = 'pearson_correlation', graph = 5, trajectories = 1000, wp = 15, dr = 'pca', core_no = 2)


# tSpace: UMAP


ts <- tSpace(df = v, K = 20, L = 15, D = 'pearson_correlation', graph = 5, trajectories = 1000, wp = 15, dr = 'umap', core_no = 2)
