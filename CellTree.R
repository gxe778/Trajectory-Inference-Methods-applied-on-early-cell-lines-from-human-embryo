#Cell Tree Gibbs: HEE

setwd("C:/Users/G177/Desktop/Trajectory Inference Repetoire")


BiocManager::install("Rgraphviz")
library(Rgraphviz)

#load HEE
load("HEE.Rda")

#load CellTree
if(!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("cellTree")

library(cellTree)

HEEA<-as.matrix(HEE@assayData$exprs)

#select K topics jusing LDA inference (to match differentiation steps)
lda.results = compute.lda(HEEA, k.topics=3:8, method="maptpx")

# alternative: using Gibbs method(much slower)
#lda.results = compute.lda(HEEA, k.topics=6, method="Gibbs")

#convert to HGNC gene names
HEEA.hgnc = HEEA

library("biomaRt") 

ensembl.ids = sapply(strsplit(rownames(HEEA), split=".",fixed=TRUE), "[", 1)

ensembl.mart = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

gene.map = getBM(attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensembl.ids, mart = ensembl.mart) 
                 

idx = match(ensembl.ids, gene.map$ensembl_gene_id) 
hgnc.ids = gene.map$hgnc_symbol[idx] 
has.hgnc.ids = !is.na(hgnc.ids)&(hgnc.ids!="") 
rownames(HEEA.hgnc)[has.hgnc.ids] = hgnc.ids[has.hgnc.ids]

HEEA_lda_model = compute.lda(HEEA.hgnc, k.topics=6)

# Number of topics of fitted model: 
print(HEEA_lda_model$K)

# Compute pairwise distance between cells based on topic distributions in the fitted model: 
dists = get.cell.dists(HEEA_lda_model)

#set phenotype/sampling data
days.factor = HEE@phenoData@data$Sample
days =levels(days.factor)[days.factor]

#compute MST:
mst.tree = compute.backbone.tree(HEEA_lda_model, days, only.mst=TRUE)

mst.tree.with.layout = ct.plot.topics(mst.tree)

# showing time point for each cell
mst.tree.with.layout = ct.plot.grouping(mst.tree)

#backbone for overall differentiation path
b.tree = compute.backbone.tree(HEEA_lda_model, days)

b.tree.with.layout = ct.plot.grouping(b.tree)

#expand beginning
b.tree = compute.backbone.tree(HEEA_lda_model, days, width.scale.factor=1.5)

b.tree.with.layout = ct.plot.grouping(b.tree)

#plot with topics assigned
b.tree.with.layout = ct.plot.topics(b.tree)

#Load GO mappings for human: 
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#GO results for each topic: Cellular Components("CC")
go.results = compute.go.enrichment(HEEA_lda_model, org.Hs.eg.db, ontology.type="CC", bonferroni.correct=TRUE, p.val.threshold=0.01)

# Print ranked table of significantly enriched terms for topic 1 # that do not appear in other topics: 
go.results$unique[[1]]


# Compute GO enrichment sets (using the Biological Process category) # for each topic and saves DAG plots to files:
go.results.bp = compute.go.enrichment(HEEA_lda_model, org.Hs.eg.db, ontology.type="BP", bonferroni.correct=TRUE, p.val.threshold=0.01, dag.file.prefix="HEEA_go")

# plot GO sub-DAG for topics 1 to 3: 
go.dag.subtree = ct.plot.go.dag(go.results, up.generations = 2, only.topics=c(1:3))

#Generate table summary of cells, ranked by tree position: 
cell.table = cell.ordering.table(b.tree)

# Print first 5 cells: 
cell.table[1:5,]


