SLICE: Determining Cell Differentiation and Lineage based on Single Cell Entropy

=================================
DESCRIPTION
=================================

Author: Minzhe Guo (minzhe.guo@cchmc.org)
Version: a12242016
License: GNU General Public License v3 <http://www.gnu.org/licenses/>.

# This is the R script of the SLICE algorithm, which utilizes single-cell RNA-seq (scRNA-seq) data to quantitatively measure cellular differentiation states based on single cell entropy and predict cell differentiation lineages via the construction of entropy directed cell trajectories. 
# To use this script, you will need the R statistical computing environment (version 3.2.2 or later) and several packages freely available through Bioconductor and CRAN, including
# * Bioconductor::Biobase, R::ggplot2, R::igraph, R::reshape2, R::entropy, R::cluster, Bioconductor::graph, Bioconductor::BioNet, R::princurve, R: lmtest, R::mgcv
# SLICE will try to resolve dependencies automatically. If the dependencies cannot be resolved, please refer to the website of each package for more information.

# SLICE is under active development. Core features have been implemented. We are improving the documentation and the visualization functions, and refining the user interfaces.
# Updates of SLICE will be distributed primarily through the SLICE website at: http://research.cchmc.org/pbge/slice.html. A step-by-step tutorial of running SLICE will also be available soon on the website. If you have any questions, please contact Dr. Minzhe Guo at minzhe.guo@cchmc.org or Dr. Yan Xu at yan.xu@cchmc.org.


# If you publish results obtained using SLICE, please cite
#   Guo M, Bao EL, Wagner M, Whitsett JA, Xu Y. 2016. SLICE: determing cell differentiation and lineage based on single cell entropy. Nucleic Acids Research. doi:10.1093/nar/gkw1278.

===================================
INSTALLATION AND RUNNING THE DEMOs
===================================

1) Download and unzip the SLICE package.
2) Open R GUI (Instructions for downloading and installation of the latest version of R computing environment can be found at http://cran.rstudio.com/).
3) Change the directory of R GUI to the directory of SLICE.
4) Run a demonstration by sourcing its R script in R GUI.
5) Analysis results are saved in the current working directory.


=================================
Files in the SLICE package
=================================

R scripts
1) slice.R: R functions to implement the SLICE algorithm.
2) AT2.R: the demonstration script of using SLICE to analyze AT2 data.
3) FB.R: the demonstration script of using SLICE to analyze FB data.
4) HEE.R: the demonstration script of using SLICE to analyze HEE data.
5) HSMM.R: the demonstration script of using SLICE to analyze HSMM data.
6) EPI.R: the demonstration script of using SLICE to analysis EPI data.

The data folder contains all the data files needed for reproducing the results in the manuscript
1) GSE52583.Rda: R data object file containing scRNA-seq data from Treutlein et al., 2014 for SLICE demonstration, including alveolar type 2 cells (AT2) from E14.5, E16.5, E18.5, and adult mouse lung, and three populations of epithelial cells at E18.5 (EPI). Expression profiles with same gene symbol annotations were averaged.
2) FB.Rda: R data object file containing scRNA-seq data of five predicted fibroblastic subtype cells (FB) from mouse lung at E16.5 (Guo et al., 2015; Du et al., 2015). Expression profiles with same gene symbol annotations were averaged.
3) HEE.Rda: R data object file containing scRNA-seq data from single cells (HEE, n=88) from seven developmental stages of human early embryo (Yan et al., 2013). Expression profiles with same gene symbol annotations were averaged.
4) HSMM.Rda: R data object file containing scRNA-seq data from differentiating human skeletal muscle myoblasts (HSMM) (Trapnell et al., 2014). Expression profiles with same gene symbol annotations were averaged.
5) mm_km.Rda: R data object file containing a pre-compiled pairwise Kappa similarity (based on GO_BP_FAT annotations downloaded from DAVID on Dec. 13, 2015) of mouse genes.  
6) hs_km.Rda: R data object file containing a pre-compiled pairwise Kappa similarity (based on GO_BP_FAT annotations downloaded from DAVID on Dec. 13, 2015) of human genes.
7) GSE52583-sig.txt: signature genes identified in the original analysis (Treutlein et al., 2014) and used in the SLICE analysis of EPI data
8) FB-sig.txt: cell type signature genes identified in the original analysis (Guo et al., 2015; Du et al., 2015) and used in the SLICE analysis of FB data
9) HSMM-markers.txt: marker genes important in myogenesis and used in the SLICE analysis of HSMM data


===================================================
Changes from the previous version a12212016
===================================================

* Updated code comments in the R scripts
* No code updates.


=================================
DEPENDENCIES
=================================

SLICE will try to resolve the dependencies automatically. 
If the dependencies cannot be resolved, please try the following scripts for installation or refer to the website of each package for more information.


if(!require(Biobase)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biobase")
}
if(!require(ggplot2)) {
    install.packages('ggplot2', dep=T)
}
if (!require(igraph)) {
    install.packages('igraph', dep=T)
}
if (!require(reshape2)) {
    install.packages('reshape2', dep=T)
}
if(!require(entropy)){
	install.packages('entropy', dep=T)
}
if(!require(graph)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("graph")
}
if(!require(BioNet)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("BioNet")
}
if(!require(cluster)){
	install.packages('cluster', dep=T)
}
if(!require(princurve)){
    install.packages('princurve', dep=T)
}
if(!require(lmtest)){
  install.packages('lmtest', dep=T)
}
if(!require(mgcv)){
  install.packages('mgcv', dep=T)
}
if(!require(gridExtra)){
  install.packages('gridExtra', dep=T)
}



require(gridExtra)
require(Biobase)
require(graph)
require(BioNet)
require(entropy)
require(cluster)
require(ggplot2)
require(grid)
require(princurve)
require(splines)
require(mgcv)
require(lmtest)
require(igraph)





