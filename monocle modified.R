#Monocle DDRTree Script



#Load Monocle:

#devtools::install_github('cole-trapnell-lab/monocle3')



library(monocle3)



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




#pre-process, partition, reduce, and cluster Monocle 3

cds <- new_cell_data_set(v)



cds <- preprocess_cds(cds, num_dim = 100)



#Reduce via UMAP

cds <- reduce_dimension(cds)



cds <- cluster_cells(cds)



cds <- learn_graph(cds)





plot_cells(cds)





partition<-colnames(v)





plot_cells(cds, group_cells_by="partition")







plot_cells(cds,
           
           color_cells_by = "partition",
           
           label_groups_by_cluster=FALSE,
           
           label_leaves=FALSE,
           
           label_branch_points=FALSE)



cds = order_cells(cds)



plot_cells(cds,
           
           color_cells_by = "pseudotime",
           
           label_cell_groups=FALSE,
           
           label_leaves=FALSE,
           
           label_branch_points=FALSE,
           
           graph_label_size=1.5)

