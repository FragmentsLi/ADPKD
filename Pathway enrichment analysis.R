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

# 2)Creating Supplymentary Fig.2
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


####################################################################
#   2. Pathway Enrichment for MNPs
####################################################################
load('MNP.RData')
#!!!只取ADPKD，15059个细胞
MNP=subset(x = MNP, subset = batch == "ADPKD")

#hallmark
vision.obj <- Vision(MNP,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=5)
vision.obj <- analyze(vision.obj)

#KEGG
vision.obj <- Vision(MNP,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=5)
vision.obj <- analyze(vision.obj)

#GOBP
vision.obj <- Vision(MNP,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)



#获得每个细胞的不同Gene Set打分
sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)

#VISION得出的score取平均值
group.by <- "subtype" 
#mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
sigScores <- expm1(x = sigScores) 
sigScores <- aggregate(sigScores, by=list(MNP@meta.data[[group.by]]), FUN="mean") 
rownames(sigScores) <- sigScores$Group.1 
sigScores <- t(sigScores[,-1]) 
head(sigScores)

#FIBROBLAST
sigScores2<-sigScores[rownames(sigScores)[grepl('FIBRO',rownames(sigScores))],]

#INTEGRIN
sigScores2<-sigScores[rownames(sigScores)[grepl('INTEGRIN',rownames(sigScores))],]

#PROLIFERATION
sigScores2<-sigScores[rownames(sigScores)[grepl('PROLIFERATION',rownames(sigScores))],]

#INFLAMMATION
sigScores2<-sigScores[rownames(sigScores)[grepl('INFLAMMA',rownames(sigScores))],]

#MIGRATION
sigScores2<-sigScores[rownames(sigScores)[grepl('MIGRATION',rownames(sigScores))],]

#ROS
sigScores2<-sigScores[rownames(sigScores)[grepl('OXIDATIVE',rownames(sigScores))],]

#DIFF
sigScores2<-sigScores[rownames(sigScores)[grepl('DIFFERENT',rownames(sigScores))],]

#DIFF
sigScores2<-sigScores[rownames(sigScores)[grepl('MACROPHAGE',rownames(sigScores))],]

#Fig：最终画图
final<-c('GOBP_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE',
         'GOBP_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE',
         'GOBP_NATURAL_KILLER_CELL_DIFFERENTIATION',
         'GOBP_B_CELL_DIFFERENTIATION',
         'GOBP_MONONUCLEAR_CELL_DIFFERENTIATION',
         'GOBP_T_CELL_DIFFERENTIATION',
         'GOBP_MAST_CELL_DIFFERENTIATION',
         'GOBP_MYOFIBROBLAST_DIFFERENTIATION',
        'GOBP_MACROPHAGE_DIFFERENTIATION',
        'GOBP_REGULATION_OF_MACROPHAGE_PROLIFERATION',
        'GOBP_REGULATION_OF_CELL_CHEMOTAXIS_TO_FIBROBLAST_GROWTH_FACTOR',
        'GOBP_FIBROBLAST_MIGRATION',
        'GOBP_MACROPHAGE_PROLIFERATION',
        'GOBP_EPITHELIAL_CELL_PROLIFERATION',
        'GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION')
sigScores2<-sigScores[final,]
rownames(sigScores2)<-c('Regulation of acute inflammatory response',
                        'Positive regulation of acute inflammatory response',
                        'Natural killer cell differentiation',
                        'B cell differentiation',
                        'Mononuclear cell differentiation',
                        'T cell differentiation',
                        'Mast cell differentiation',
                        'Myofibroblast differentiation',
                        'Macrophage_differentiation',
                        'Regulation of macrophage proliferation',
                        'Regulation of cell chemotaxis to fibroblast growth factor',
                        'Fibroblast migration',
                        'Macrophage proliferation',
                        'Epithelial cell proliferation',
                        'Positive regulation of epithelial cell proliferation')

#   Fig: GO结果热图
p <- pheatmap::pheatmap(sigScores2,scale = "row",
                        cluster_cols = F,
                        cluster_rows = F,
                        fontsize=8,
                        border_color = "grey",
                        cellwidth = 12, 
                        cellheight = 12,
                        angle_col=90,
                        )


png(filename = "ADPKD_MNP_VISION_GOBP_final.png",width = 120, height = 120, 
    units = "mm", res = 300)
p
dev.off()

#write.csv(sigScores,'VISION_heatmap_KEGG.csv',quote=F)

#   Fig:hallmark结果热图
p <- pheatmap::pheatmap(sigScores,scale = "row",cluster_cols = F,
                        fontsize=8,
                        border_color = "grey",
                        cellwidth = 12, 
                        cellheight = 12,
                        annotation_colors=ann_colors,angle_col=90)

png(filename = "ADPKD_MNP_VISION_hallmark_all.png",width = 200, height = 290, 
    units = "mm", res = 300)
p
dev.off()

