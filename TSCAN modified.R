#TSCAN script

#Load Data
devtools::install_github("TSCAN")
library(TSCAN)

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


#load/process data
data<-v

data <- data[rowSums(data) > 0,]

clures <- hclust(dist(data))

cluster <- cutree(clures,0.05*nrow(data))

aggdata <- aggregate(data,list(cluster),mean)

aggdata <- aggdata[,-1]


#kmeans

HSMMkmeans <- kmeans(aggdata, centers=7)



HSMMkmeansorder <- TSCANorder(HSMMkmeansorder)

#waterfall example
HSMMwaterfallorder <- pseudotimeprog.foo(aggdata,color="black") 

