#HEE: FORKS

#load HEE
load("HEE.Rda")

#install mfa
vignette('introduction_to_mfa')

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("mfa")


#load mfa
library(mfa)

#create data
HEEA<-as.matrix(HEE@assayData$exprs)

HEEA<-t(HEE@assayData$exprs)

#remove zero variance
HEEA<-HEEA[ , apply(HEEA, 2, var) != 0]


#run mfa
m <- mfa(HEEA)

print(m)

plot_mfa_trace(m)

ms <- summary(m)
print(head(ms))

qplot(rownames(HEEA), ms$pseudotime, color = factor(rownames(HEEA))) +
  xlab('True pseudotime') + ylab('Inferred pseudotime') +
  scale_color_discrete(name = 'True\nbranch')
