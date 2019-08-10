#SCORPIOUS script

#https://github.com/rcannood/SCORPIUS


devtools::install_github("rcannood/SCORPIUS", build_vignettes = TRUE)

library(SCORPIUS)
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


#create groupings
groups<- c(colnames(v))

#Dim. reduction
space <- reduce_dimensionality(v, "spearman")
draw_trajectory_plot(space,contour = TRUE)

#Inference
traj <- infer_trajectory(space)
draw_trajectory_plot(space)

draw_trajectory_plot(traj$path, contour = TRUE)

gimp <- gene_importances(v,
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000)

gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(expression[,gene_sel])

modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules)


