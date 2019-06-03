#FATEID: HEE




#load HEE
load("HEE.Rda")

#load FATEID
install.packages("devtools")
library(devtools)
install_github("dgrun/FateID")
library(FateID)

x<-HEE$Sample
y<-HEE$Cell
x<-as.numeric(x)
y<-as.numeric(y)
x<-as.matrix(x)
y<-as.matrix(y)

#determine clusters
#end stages of cell development

tar<-HEE@assayData$exprs


#Fate bias

fb  <- fateBias(x, y, tar, z=NULL, minnr=5, minnrh=10, adapt=TRUE, confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)


dr  <- compdr(x, z=NULL, m=c("tsne","cmd","dm","lle","umap"), k=c(2,3), lle.n=30, dm.distance="euclidean", tsne.perplexity=20, seed=12345)

plotFateMap(y,dr,k=2,m="tsne")

g=rownames(x)

pr <-plotFateMap(y, dr, k=2, m="tsne", x=g, prc=TRUE)

pr <-plotFateMap(y, dr, k=2, m="tsne", x=g, trthr=.33, prc=TRUE)


pr  <- prcurve(y,dr,k=2,m="tsne",trthr=0.4,start=3)
