load('adpkd_all.RData')
library(Seurat)
library(SingleR)
library(celldex)

#   Annotation
Idents(adpkd_qc) <- 'RNA_snn_res.0.7'
type<-c('PT','Tcell','MNP','MNP','Tcell',#0,1,2,3,4
        'FIB','Bcell','Tcell','PT','FIB',#5,6,7,8,9
        'PT','MAST','EPI','FIB','Tcell',#10,11,12,13,14
        'ENDO','EPI','EPI','Bcell','MNP',#15,16,17,18,19
        'Tcell','MNP','DT','FIB','EPI',#20,21,22,23,24
        'MAST','Bcell')#25,26

names(type) <- levels(adpkd_qc) 
adpkd_qc <- RenameIdents(adpkd_qc, type)       #修改Idents
adpkd_qc$celltype <- Idents(adpkd_qc)
