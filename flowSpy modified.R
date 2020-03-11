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


   save(TID,file="TID.Rda")



   write.xlsx2(TID, file="TID.xlsx")




counts<-as.matrix(TID)



v<-log2(counts+1)


#save: excel and R:


   save(v,file="v.Rda")


   write.xlsx2(v, file="v.xlsx")



#create a SingleCellExperiment object sce


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))


#preprocessing

meta.data <- data.frame(cell = rownames(vv), stage =rownames(vv))

meta.data$stage <- factor(as.character(meta.data$stage), levels = rownames(vv))

markers<-colnames(vv)


fspy <- createFSPY(raw.data = vv, markers = markers,
                   meta.data = meta.data,
                   normalization.metho = "log2",
                   verbose = TRUE)


fspy <- runCluster(fspy, cluster.method = "som")

fspy <- processingCluster(fspy)

fspy <- runFastPCA(fspy)

#run t-Distributed Stochastic Neighbor Embedding (tSNE)

fspy <- runTSNE(fspy, perplexity = 20)

fspy <- runUMAP(fspy)

fspy <- buildTree(fspy, dim.type = "tsne", dim.use = 1:2)


#define trejectory root and stem


diff.list <- runDiff(fspy)

fspy <- defRootCells(fspy, root.cells = c("Oocyte #1(RPKM)"))

fspy <- defLeafCells(fspy, leaf.cells = c("Late blastocyst #3 -Cell#8(RPKM)"),verbose = TRUE)


#run pseudotime

fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "raw")

plotPseudotimeTraj(fspy, var.cols = FALSE, markers=markers) 

plotPseudotimeTraj(fspy, cutoff = .05, var.cols = TRUE, markers=markers)


#Plot 2D tSNE; cells colored by cluster id

plot2D(fspy, item.use = c("tSNE_1", "tSNE_2"), color.by = "cluster.id", 
       
       alpha = 1, main = "tSNE", category = "categorical", show.cluser.id = TRUE)

#plot 2D UMAP; cells colored by cluster id

plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id", 
       
       alpha = 1, main = "UMAP", category = "categorical", show.cluser.id = TRUE)




#UMAP plot colored by pseudotime

plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), category = "numeric",
       
       size = 1, color.by = "pseudotime") + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))



#tSNE plot colored by pseudotime

plot2D(fspy, item.use = c("tSNE_1", "tSNE_2"), category = "numeric",
       
       size = 1, color.by = "pseudotime") + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))



#trajectory value

plotPseudotimeTraj(fspy, var.cols = TRUE) + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))


plotPseudotimeTraj(fspy, cutoff = 0.05, var.cols = TRUE) + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))


