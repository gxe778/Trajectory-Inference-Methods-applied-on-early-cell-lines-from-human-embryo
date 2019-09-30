#flowSpy TI


#Install repository and load flowSpy
devtools::install_github("JhuangLab/flowSpy")

library(flowSpy)


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

vv<-t(v)

sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))


#Preprocessing

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

# run t-Distributed Stochastic Neighbor Embedding (tSNE)

fspy <- runTSNE(fspy, perplexity = 20)

fspy <- runUMAP(fspy)

fspy <- buildTree(fspy, dim.type = "tsne", dim.use = 1:2)


#Define Traj. root and stem


diff.list <- runDiff(fspy)

fspy <- defRootCells(fspy, root.cells = c("Oocyte #1(RPKM)"))

fspy <- defLeafCells(fspy, leaf.cells = c("Late blastocyst #3 -Cell#8(RPKM)"),verbose = TRUE)


# run pseudotime

fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "raw")

plotPseudotimeTraj(fspy, var.cols = FALSE, markers=markers) 

plotPseudotimeTraj(fspy, cutoff = .05, var.cols = TRUE, markers=markers)


# Plot 2D tSNE. And cells are colored by cluster id

plot2D(fspy, item.use = c("tSNE_1", "tSNE_2"), color.by = "cluster.id", 
       
       alpha = 1, main = "tSNE", category = "categorical", show.cluser.id = TRUE)

# Plot 2D UMAP. And cells are colored by cluster id

plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id", 
       
       alpha = 1, main = "UMAP", category = "categorical", show.cluser.id = TRUE)


# Tree plot

plotTree(fspy, color.by = "D0.percent", show.node.name = TRUE, cex.size = 1) + 
  
  scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))


#Cluster plot

plotCluster(fspy, item.use = c("tSNE_1", "tSNE_2"), category = "numeric",
            
            size = 100, color.by = "CD45RA") + 
  
  scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))


# UMAP plot colored by pseudotime

plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), category = "numeric",
       
       size = 1, color.by = "pseudotime") + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))



# tSNE plot colored by pseudotime

plot2D(fspy, item.use = c("tSNE_1", "tSNE_2"), category = "numeric",
       
       size = 1, color.by = "pseudotime") + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))



# trajectory value

plotPseudotimeTraj(fspy, var.cols = TRUE) + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))


plotPseudotimeTraj(fspy, cutoff = 0.05, var.cols = TRUE) + 
  
  scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))


