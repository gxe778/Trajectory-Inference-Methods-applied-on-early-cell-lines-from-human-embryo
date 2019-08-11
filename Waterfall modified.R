#Waterfall Script

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


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))





all <-v


# Focusing on entire frame

c <- cor(all, method="pearson")
d <- dist(c)
hr <- hclust(d, method = "ward", members=NULL)
plot(hr, hang = 0.1,labels=colnames(all))
plot(hr, hang = 0.1,labels=cutree(hr,k=7)) 

pca18 <- prcomp(as.data.frame(t(all)), cor=T)

gene_list <-unique(c(head(names(sort(pca18$rotation[,1])),18),tail(names(sort(pca18$rotation[,1])),18),head(names(sort(pca18$rotation[,2])),18),tail(names(sort(pca18$rotation[,2])),18),head(names(sort(pca18$rotation[,3])),18),tail(names(sort(pca18$rotation[,3])),18),head(names(sort(pca18$rotation[,4])),18),tail(names(sort(pca18$rotation[,4])),18)))

all <-t(scale(t(all[match(gene_list,rownames(all)),])))
pca18 <- prcomp(as.data.frame(t(all)), cor=T)

branch <-cutree(hr,k=5)
all <-rep("#BFBFBF30",ncol(all))
all[which(branch==1)] <-"#4882C3" 
all[which(branch==2)] <-"#13751B" 
all[which(branch==3)] <-"#F26A6A" 
all[which(branch==4)] <-"#FF6A00" 
all[which(branch==5)] <-"#3AC4D6" 
all[which(branch==6)] <-"##0000FF"
all[which(branch==7)] <-"#FF0000" 
names(all) <-colnames(all)

plot(pca18$x[,1:90],col=all,pch=19) # Figure M3B


# 3. O:B
gene_name=as.matrix(all)
gene_name=gene_name
plot(pca18$x[,1:90])

# 4. O:Z
gene_name=as.matrix(all)
plot(pca18$x[,1:2])

###Voronoi plot ###
all <-v

c <- cor(as.matrix(all), method="pearson")
d <- dist(c)
hr <- hclust(d, method = "ward", members=NULL)
plot(hr, hang = 0.1)
d <-dist(t(all))
plot(hr,labels=cutree(hr,k=5))
cutree(hr,h=10)

all<-as.matrix(TID)
d<-dist(t(all))
r <- sammon(d,niter=500000,tol=1e-500)
x <- r$points
library(deldir)
sammon.dd <-deldir(x[,1],x[,2])
plot(tile.list(sammon.dd),showpoints=FALSE)

r <- sammon(d,niter=50000000,tol=1e-500)
x <- r$points
library(deldir)
sammon.dd <-deldir(x[,1],x[,2])
plot(tile.list(sammon.dd),showpoints=FALSE)


pca18 <- prcomp(as.data.frame(t(all)))
plot(pca18$x[,1:90],pch=19)
text(pca18$x[,1:90],labels=branch <-branch[match(colnames(all),names(branch))])
text(pca18$x[,1:90],labels=colnames(all))

all.col <-colnames(all)

# Visualize Z and B routes
pca18 <- prcomp(as.data.frame(t(all[,1:90])))
plot(pca18,pch=19)
text(pca18$x[,1:90],labels=names(branch))

###Voronoi plot ###
c <- cor(all, method="pearson")
d <- dist(c)
hr <- hclust(d, method = "ward", members=NULL)
plot(hr, hang = 0.1); 
plot(hr,labels=colnames(all)) #;plot(hr,labels=cutree(hr,h=1))
d <-dist(t(all))
plot(hr,labels=cutree(hr,h=5))
cutree(hr,h=10)
r <- sammon(d,niter=500000,tol=1e-500)
x <- r$points
library(deldir)
sammon.dd <-deldir(x[,1],x[,2])
plot(tile.list(sammon.dd),showpoints=FALSE)




