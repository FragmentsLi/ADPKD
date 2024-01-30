library(Seurat)
load('adpkd_all.RData')
source('.../DMC-master/R/create_sig_matrix.R')

# Get normalized matrix
df <- as.data.frame(GetAssayData(adpkd_qc, assay = "RNA", slot = "counts"))
# Retain HVGs
df<-df[VariableFeatures(object = adpkd_qc),]

pheno<-data.frame(cell_type=Idents(adpkd_qc))

adpkd_qc=0

# Create signature matrix
options(mc.cores=10)
sig_matrix<-create_sig_matrix(df,pheno)

