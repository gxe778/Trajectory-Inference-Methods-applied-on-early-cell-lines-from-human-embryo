#CellTrails Script


#Load Package
BiocManager::install("CellTrails")


library(CellTrails)


#Import excel file and manipulate
TI_Data_Yan_et_al <- read_excel("TI Data Yan et al.xlsx")


TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)


TID<-TI_Data_Yan_et_al


rownames(TID)<-TI_Data_Yan_et_al[,1]


TID<-TID[,-c(1:2)]


#save data in excel and R
save(TID, file="TID.xlsx")

save(TID,file="TID.Rda")


#Create SingleCell Object
counts<-as.matrix(TID)

v<-log2(counts+1)


# create a SingleCellExperiment object sce
sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))
                           
#select features of Trajectory
showTrajInfo(sce)

tfeat <- filterTrajFeaturesByDL(sce, threshold=NULL)


tfeat


trajFeatureNames(sce) <- tfeat



showTrajInfo(sce)






#Dimension Reduction
se <- embedSamples(sce)

names(se)

d <- findSpectrum(se$eigenvalues, frac=100)

latentSpace(sce) <- se$components[,d]

showTrajInfo(sce)

featureNames(sce)

plotManifold(sce, color_by="featurename", name="ACAP2")

cl <- findStates(sce, min_size=0.01, min_feat=5, max_pval=1e-4, min_fc=2)


head(cl)



states(sce) <- cl

sce <- connectStates(sce, l=10)

trajComponents(sce)

plotStateTrajectory(sce, color_by="featureName", name="ACAP2", 
                    component=1, point_size=1.5, label_offset=4)


gp<-plotStateTrajectory(sce, color_by="phenoName", name="state", 
                    component=1, point_size=1.5, label_offset=4)


stateTrajLayout(sce) <- gp


plotStateTrajectory(sce, color_by="featureName", name="ACAP2", 
                    component=1, point_size=5)

sce <- selectTrajectory(sce, component=1)

sce_subset <- sce[, trajSampleNames(sce)]

plotStateSize(sce)

sce <- fitTrajectory(sce)

showTrajInfo(sce)

sce <- fitTrajectory(sce)

sce <- fitTrajectory(sce)

plotTrajectoryFit(sce)
