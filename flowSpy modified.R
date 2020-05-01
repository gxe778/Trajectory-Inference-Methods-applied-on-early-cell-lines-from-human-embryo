                                  #flowSpy 


#devtools::install_github("JhuangLab/flowSpy")

library(flowSpy)


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


#save: excel and R:


                      #save(v,file="v.Rda")


                      #write.xlsx2(v, file="v.xlsx")



#create a SingleCellExperiment object sce


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))


#preprocessing

meta.data <- data.frame(cell = rownames(v), stage = rownames(v))


meta.data$stage <- factor(as.character(meta.data$stage), levels = rownames(v))

markers<-colnames(v)


fspy <- createFSPY(raw.data = v, markers = markers,
                   meta.data = meta.data,
                   normalization.metho = "log2",
                   verbose = TRUE)


fspy <- runCluster(fspy, cluster.method = "som")

fspy <- processingCluster(fspy)


fspy <- runUMAP(fspy)

fspy <- buildTree(fspy, dim.use = 1:2)




diff.list <- runDiff(fspy)



#run pseudotime

fspy <- runPseudotime(fspy, verbose = TRUE, mode= "undirected", dim.use = 1:2, dim.type = "umap")


plotPseudotimeTraj(fspy, var.cols = FALSE, markers=markers) 





#plot 2D UMAP: cells colored by cluster id


plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id", alpha = 1, main = "UMAP", category = "categorical", show.cluser.id = TRUE)




#UMAP plot colored by pseudotime


plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), category = "numeric", size = 1, color.by = "pseudotime", scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0")))






#trajectory value



plotPseudotimeTraj(fspy, var.cols = TRUE) 


