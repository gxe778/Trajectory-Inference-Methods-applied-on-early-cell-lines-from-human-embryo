


#Calista Script



#devtools::install_github("CABSEL/CALISTA/CALISTA-R/calista")





#Load Data

library(calista)





#Load Data

TI_Data_Yan_et_al <- read_excel("TI_Yan.xlsx")



TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)





TID<-TI_Data_Yan_et_al





rownames(TID)<-TI_Data_Yan_et_al[,1]





TID<-TID[,-c(1:2)]





# save data in excel and R:


# save(TID,file="TID.Rda")



# write.xlsx2(TID, file="TID.xlsx")




counts<-as.matrix(TID)



v<-log2(counts+1)


# save: excel and R:


# save(v,file="v.Rda")


# write.xlsx2(v, file="v.xlsx")



# create a SingleCellExperiment object sce


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))





#Specify data types and settings for pre-processing**

# 1- for scRT-qPCR  (CT values)

# 2- for scRT-qPCR  (Expression values - Num. mRNA molecules)

# 3- for scRNA-seq  (Expression values - e.g log(TPM+1) or log(RPKM+1))

# 4- for scRNA-seq  (Expression values - e.g TPM+1 or RPKM)



INPUTS=list()

INPUTS$data_type=3



# 1- Rows= cells and Columns= genes with time/stage info in the last column  

# 2- Rows= genes and Columns= cells with time/stage info in the first row

# 3- Rows= cells and Columns= genes without time info

# 4- Rows= genes and Columns= cells without time info

# 5- Manual selection from data table

INPUTS$format_data=4



# When cells come from different capture times or stages, users can select 

# cells from specific time or stage for further analysis by CALISTA. 

# For instance, considering a dataset with cells taken from 4 time points, 

# snapshots can be selected, as follows:

# INPUTS$data_selection = integer() or c(1:4) for all time points/cells 

# INPUTS$data_selection = 1          for cells the first time point

# INPUTS$data_selection = c(1 3 4)    for cells 1st, 3rd and 4th time points

INPUTS$data_selection=integer() 



# Users can exclude genes with a certain percentage of zeros. 

# INPUTS$perczeros_genes = 100 (recommended)

INPUTS$perczeros_genes=100



# Remove cells with 100% of zeros

# Users can exclude cells with more than a certain percentage of zeros.

# INPUTS$perczeros_cells = 100 (recommended)



INPUTS$perczeros_cells=100



# Users can exclude cells from further analysis. 

# 1- Remove cells based on their indices in the expression matrix. 

#    Indices need to be uploaded as a separate csv file. 

# 0- No cells deletion



INPUTS$cells_2_cut=0 



# Retain only top X the most variable genes with X=min(200, INPUTS.perc_top_genes * num of cells/100, num of genes)

INPUTS$perc_top_genes=10 



# Specify single-cell clustering settings

# 1- select the number of clusters by eigengap plot

# 0- define the number of clusters 



INPUTS$optimize=1 



# 1- Use parallel processing (number of cores available - 1)  

# 0- Do not use processing



INPUTS$parallel=1



# Number of clustering runs (for greedy algorithm)

# INPUTS$runs = 50; 

INPUTS$runs=50 



# Maximum number of iterations in the greedy algorithm 

# INPUTS$max_iter = 100; 



INPUTS$max_iter=100 



# 'hierarchical'- hierachical clustering of consensus matrix

# 'kmedoids'-  kmedoids for the clustering of consensus matrix

INPUTS$Cluster='kmedoids' 



#Set threshold for transition genes determination to 50%

INPUTS$thr_transition_genes=75 



#v, in this script, has been imported into excel and saved as an xlsv file. Load "v.xlsv"(listed in repository)

calista::calista(INPUTS)



