####################################################################
#   1. Pathway Enrichment for Major Cell Lineages
####################################################################
library(VISION)
library(Seurat)
load('adpkd_all.RData')
# 1)VISION
vision.obj <- Vision(adpkd_qc,
              signatures = c(".../msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
vision.obj=0

# 2)Mean value
group.by <- "celltype2" 
sigScores <- expm1(x = sigScores) 
sigScores <- aggregate(sigScores, by=list(adpkd_qc@meta.data[[group.by]]), FUN="mean") 
rownames(sigScores) <- sigScores$Group.1 
sigScores <- t(sigScores[,-1]) 
head(sigScores)

adpkd<-c('PT_ADPKD','Tcell_ADPKD','MNP_ADPKD',
        'FIB_ADPKD','Bcell_ADPKD','MAST_ADPKD','EPI_ADPKD',
        'ENDO_ADPKD','DT_ADPKD')
control<-c('PT_Normal','Tcell_Normal','MNP_Normal',
        'FIB_Normal','Bcell_Normal','MAST_Normal','EPI_Normal',
        'ENDO_Normal','DT_Normal')
sigScores<-sigScores[,c(adpkd,control)]
