#FATE ID SCRIPT
#https://github.com/dgrun/FateID

install_github("dgrun/FateID")

library(FateID)

#load Data
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


#FATE Analysis
tar<-c(v[, c(10,30,50)])

tar <- tar

x <- v

y<-v

fb  <- fateBias(x,y, tar, z=NULL, minnr=5, minnrh=10, adapt=TRUE, confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)

dr  <- compdr(x, z=NULL, m=c("tsne","cmd","dm","lle","umap"), k=c(2,3), lle.n=30, dm.sigma="local", dm.distance="euclidean", tsne.perplexity=30, seed=12345)