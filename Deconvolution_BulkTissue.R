## RNA-Seq deconvolution ---
### 

# environment
# library("purrr")
library("tibble")
library("dplyr")
library("tidyr")
library("tidyverse")
select <- dplyr::select

source("Deconvolution_functions.R")

writeTXTfile <- function(variable, filename){
  write.table(cbind(rownames(variable),variable), file=paste0(filename), row.names=FALSE, quote=FALSE, sep="\t")
}
readTXTfile = function(filename,rownames){
  if(rownames==0){read.csv(file=filename, sep="\t", header=T, stringsAsFactors = F)}else{
    read.csv(file=filename, sep="\t", header=T, stringsAsFactors = F, row.names=rownames)
  }}


basepath = "/crex/proj/sctatlas/nobackup/private/mengnanCopy/scRNA_atlas_v4/dwls/"
# setwd(basepath)

# basepath=file.path(input.path,"rna_single_cell_read_count")
dataBulk = readTXTfile(file.path(basepath,"data","hpa_dataBulk_hpaonly.txt"),1)
celltype_enriched_genes = readTXTfile(file.path(basepath,"data","gene_classification_long.txt"),1)
celltype_enriched_genes = celltype_enriched_genes %>% 
  filter(specificity_category=="Cell type enriched",dataset=="singlecell") %>% 
  distinct(Ensembl) %>% 
  pull(Ensembl)

  # read single cell data
  dataSC = readTXTfile(file.path(basepath,"data","hpa_sc_celltype_ntpm.txt"),1)
  dataSC = dataSC %>% column_to_rownames("Gene")
  # deconvolution
  # load("/Users/mengnan.shi/Documents/project/single_cell_atlas/v4/results/deconvolution/Sig.RData")
  # Signature <- buildSignatureMatrixMAST(dataSC, labels, "results_hpa_con")
  Signature = as.matrix(dataSC[celltype_enriched_genes,])
  
  allCounts_DWLS<-NULL
  allCounts_OLS<-NULL
  allCounts_SVR<-NULL
  for(j in 1:(dim(dataBulk)[2])){
    S<-log2(Signature+1)
    Bulk<-dataBulk[,j] # extract every sample
    names(Bulk)<-rownames(dataBulk)
    Genes<-intersect(rownames(S),names(Bulk))
    B<-log2(Bulk[Genes]+1)
    S<-S[Genes,]
    solOLS<-solveOLS(S,B)
    solDWLS<-solveDampenedWLS(S,B)
    solSVR<-solveSVR(S,B)
    
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
    allCounts_OLS<-cbind(allCounts_OLS,solOLS)
    allCounts_SVR<-cbind(allCounts_SVR,solSVR)
  }
  
  colnames(allCounts_DWLS) <- colnames(dataBulk)
  colnames(allCounts_OLS) <- colnames(dataBulk)
  colnames(allCounts_SVR) <- colnames(dataBulk)
  
  writeTXTfile(allCounts_DWLS, paste0(basepath,"results_hpa/","allCounts_DWLS.txt"))
  writeTXTfile(allCounts_OLS, paste0(basepath,"results_hpa/","allCounts_OLS.txt"))
  writeTXTfile(allCounts_SVR, paste0(basepath,"results_hpa/","allCounts_SVR.txt"))
  


