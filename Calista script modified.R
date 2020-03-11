                                #Calista 



#devtools::install_github("CABSEL/CALISTA/CALISTA-R/calista")

library(calista)


#load data

TI_Data_Yan_et_al <- read_excel("TI_Yan.xlsx")



TI_Data_Yan_et_al<-as.data.frame(TI_Data_Yan_et_al, stringsAsFactors = FALSE)





TID<-TI_Data_Yan_et_al





rownames(TID)<-TI_Data_Yan_et_al[,1]





TID<-TID[,-c(1:2)]





#save data in excel and R:


                      #save(TID,file="TID.Rda")



                      #write.xlsx2(TID, file="TID.xlsx")




counts<-as.matrix(TID)



v<-log2(counts+1)


#save: excel and R:


                            #save(v,file="v.Rda")


                            #write.xlsx2(v, file="v.xlsx")



#create a SingleCellExperiment object sce


sce <- SingleCellExperiment(assay = list(counts = counts, logcounts=v))





#specify data types and settings for pre-processing#

#1- scRT-qPCR  (CT values)

#2- scRT-qPCR  (Expression values - Num. mRNA molecules)

#3- scRNA-seq  (Expression values - e.g log(TPM+1) or log(RPKM+1))

#4- scRNA-seq  (Expression values - e.g TPM+1 or RPKM)



INPUTS=list()

INPUTS$data_type=3



#1- rows= cells and columns= genes with time/stage info in the last column  

#2- rows= genes and columns= cells with time/stage info in the first row

#3- rows= cells and columns= genes without time info

#4- rows= genes and columns= cells without time info

#5- manual selection from data table

INPUTS$format_data=4



#select cells from specific time/stage: example: cells from 4 time points:


                                  #all time points/cells:

                                       #INPUTS$data_selection = integer() or c(1:4)

                                  #cells the first time point:

                                       #INPUTS$data_selection = 1  

                                  #1st, 3rd and 4th time points:

                                        #INPUTS$data_selection = c(1 3 4)     

INPUTS$data_selection=integer() 



#exclude genes with certain percentage of zeros

INPUTS$perczeros_genes=100



#remove cells with 100% of zeros:


INPUTS$perczeros_cells=100


#filter:

           #1- remove cells based on indices in expression matrix; indices uploaded as separate csv file. 

            #0- no cells deletion



INPUTS$cells_2_cut=0 



#retain only top X variable genes: X=min(200, INPUTS.perc_top_genes * num of cells/100, num of genes)

INPUTS$perc_top_genes=10 



#specify single-cell clustering settings:

             #1- select the number of clusters by eigengap plot

             #0- define the number of clusters 



INPUTS$optimize=1 


#parallel processing

               #1- parallel processing (number of cores available - 1)  

               #0- do not use processing



INPUTS$parallel=1



#number of clustering runs (greedy algorithm)


INPUTS$runs=50 



# Maximum number of iterations in greedy algorithm 


INPUTS$max_iter=100 



#'hierarchical'- hierachical clustering of consensus matrix

#'kmedoids'-  kmedoids clustering of consensus matrix

INPUTS$Cluster='kmedoids' 



#threshold for transition genes

INPUTS$thr_transition_genes=75 



#v.xlsx: excel file of expression matrix for calista.

calista::calista(INPUTS)



