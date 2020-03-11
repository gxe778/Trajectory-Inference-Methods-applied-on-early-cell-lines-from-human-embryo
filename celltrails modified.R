                                              #CellTrails 


BiocManager::install("CellTrails")

library(CellTrails)



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




#select features of trajectory



showTrajInfo(sce)





tfeat <- filterTrajFeaturesByDL(sce, threshold=NULL)




tfeat




trajFeatureNames(sce) <- tfeat




showTrajInfo(sce)





#dimension reduction



se <- embedSamples(sce)




names(se)




d <- findSpectrum(se$eigenvalues, frac=100)




latentSpace(sce) <- se$components[,d]



showTrajInfo(sce)




featureNames(sce)




plotManifold(sce, color_by="featurename", name="ACAP2")




trajSampleNames(sce)





states(sce) <- trajSampleNames(sce)




sce <- connectStates(sce, l=10)




trajComponents(sce)




showTrajInfo(sce)





#use first component

gp<-plotStateTrajectory(sce, color_by="phenoName", name="state", component=1,
                        
                        
                        
                        point_size=1.5, label_offset=4)





#use last component



gp2<-plotStateTrajectory(sce, color_by="phenoName", name="state", 
                         
                         
                         
                         component=5, point_size=.5, label_offset=4)





gp



gp2



scec <- selectTrajectory(sce, component=1)



scec2 <- selectTrajectory(sce, component=5)





plotStateSize(scec)

plotStateSize(scec2)






scec <- fitTrajectory(scec)

scec2 <- fitTrajectory(scec2)
