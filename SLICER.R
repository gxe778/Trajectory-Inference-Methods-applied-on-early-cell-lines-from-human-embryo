                                                      #SLICER Script

install.packages("SLICER")

library(SLICER)


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


#SLICER: First 300 genes

genes = rownames(v)[1:300]
k = select_k(v[1:300,], kmin=5)
traj_lle<- lle::lle(v[1:300,], m=7, k)$X
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = 1
cells_ordered = cell_order(traj_graph, start)
branches = assign_branches(traj_graph,start, min_branch_len = 10)

dists<-process_distance(traj_graph, start)

graph_process_distance(traj_graph,traj_lle,start)

graph_gene(v,traj_lle,colnames(v),1)

