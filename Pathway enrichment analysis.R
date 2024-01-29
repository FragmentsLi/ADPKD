library(VISION)
library(Seurat)

####################################################################
#   1. Pathway Enrichment for Major Cell Lineages
####################################################################
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

# 2)Creating Supplementary Fig.2
annotation_col<-data.frame(Group=factor(c(rep('ADPKD',9),rep('Normal',9))),
                           Type=factor(rep(c('PT','Tcell','MNP',
        'FIB','Bcell','MAST','EPI',
        'ENDO','DT'),2)))
rownames(annotation_col)<-colnames(sigScores)

library(ggsci)
ann_colors<-list(Group=c('ADPKD'='#f47720','Normal'='#459943'),
                       Type=c('PT'=pal_npg()(9)[1],'Tcell'=pal_npg()(9)[2],
                             'MNP'=pal_npg()(9)[3],'FIB'=pal_npg()(9)[4],
                             'Bcell'=pal_npg()(9)[5],'MAST'=pal_npg()(9)[6],
                             'EPI'=pal_npg()(9)[7],'ENDO'=pal_npg()(9)[8],
                             'DT'=pal_npg()(9)[9]))

p <- pheatmap::pheatmap(sigScores,scale = "row",cluster_cols = F,
                        annotation_col=annotation_col,
                        fontsize=8,
                        border_color = "grey",
                        cellwidth = 12, 
                        cellheight = 12,
                        annotation_colors=ann_colors,angle_col=90,
                        labels_col=rep(c('PT','Tcell','MNP',
        'FIB','Bcell','MAST','EPI',
        'ENDO','DT'),2))

png(filename = "ADPKD_VISION_hallmark_all.png",width = 230, height = 280, 
    units = "mm", res = 300)
p
dev.off()

