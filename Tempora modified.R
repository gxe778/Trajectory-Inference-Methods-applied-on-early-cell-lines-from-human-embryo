                                                  #Tempora

devtools::install_github("BaderLab/Tempora")

library(Tempora)


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



#create Seurat object 


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))

S<- CreateSeuratObject(counts = counts, project = "YAN", min.cells = 3, min.features = 200)


dense.size <- object.size(S)

sparse.size <- object.size(S)

dense.size/sparse.size

S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-")

percent.mt<-S[["percent.mt"]]

head(S, 5)

VlnPlot(S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



plot1 <- FeatureScatter(S, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

S[["RNA"]]@data



S <- SubsetData(S)

S<- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)





pbmc<-S

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


pbmc <- ScaleData(pbmc)

pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")




pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)




VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")


DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)






pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)




JackStrawPlot(pbmc, dims = 1:15)


ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)



#cluster IDs of the first 5 cells
head(Idents(pbmc), 5)



DimPlot(pbmc, reduction = "pca", label = TRUE)





#meta.data

pbmc@meta.data$Timepoints<-rownames(pbmc@meta.data)
Timepoints<-pbmc@meta.data$Timepoints

timepoint_order<-Timepoints

cluster_labels <- pbmc@meta.data$seurat_clusters
cluster_labels<-cluster_labels

pbmc@meta.data$Clusters<-pbmc@meta.data$RNA_snn_res.0.5
Clusters <- pbmc@meta.data$RNA_snn_res.0.5

celltype_markers<-NULL


#Tempora object 

T <- ImportSeuratObject(pbmc, "RNA", clusters = Clusters,
                                     timepoints = Timepoints, 
                                     cluster_labels = levels(cluster_labels),
                                     timepoint_order = timepoint_order)



#download Human Enrichment Set: Human_GO_AllPathways_no_GO_iea_March_01_2020_symbol.gmt:   

                     #http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
  




#estimate pathway enrichment profiles of clusters; Human Enrichment Set within working directory

cortex_tempora <- CalculatePWProfiles(T, 
                                      gmt_path =".....",
                                      method="gsva", min.sz = 5, max.sz = 200, parallel.sz = 1)

#build trajectory 

cortex_tempora <- BuildTrajectory(T, n_pcs = ...., difference_threshold = 0.01)


#visualize trajectory

cortex_tempora <- PlotTrajectory(T)

#fit GAMs on pathway enrichment profile

cortex_tempora <- IdentifyVaryingPWs(T, pval_threshold = 0.05)
