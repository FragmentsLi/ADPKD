library(VISION)
library(Seurat)
load('adpkd_all.RData')


#运行VISION，！注意很耗内存，耗时，10个核约花2h+
#1.做Hallmark的富集分析
vision.obj <- Vision(adpkd_qc,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

#2.做KEGG的富集分析
vision.obj <- Vision(adpkd_qc,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

#3.做PID的富集分析
vision.obj <- Vision(adpkd_qc,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c2.cp.pid.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

#4.做GO的富集分析(BP太大暂时不做)
vision.obj <- Vision(adpkd_qc,
              signatures = c("/data4/JN/PRJNA679848/VISION/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c5.go.mf.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)


#获得每个细胞的不同Gene Set打分
sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
#清除释放内存
vision.obj=0

#VISION得出的score取平均值
group.by <- "celltype2" 
#mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
sigScores <- expm1(x = sigScores) 
sigScores <- aggregate(sigScores, by=list(adpkd_qc@meta.data[[group.by]]), FUN="mean") 
rownames(sigScores) <- sigScores$Group.1 
sigScores <- t(sigScores[,-1]) 
head(sigScores)

write.csv(sigScores,'VISION_heatmap_hallmark.csv',quote=F)

#给细胞类型列序
sigScores<-read.csv('VISION_heatmap_hallmark.csv',header=T,check.names=F,row.names=1)
adpkd<-c('PT_ADPKD','Tcell_ADPKD','MNP_ADPKD',
        'FIB_ADPKD','Bcell_ADPKD','MAST_ADPKD','EPI_ADPKD',
        'ENDO_ADPKD','DT_ADPKD')
control<-c('PT_Normal','Tcell_Normal','MNP_Normal',
        'FIB_Normal','Bcell_Normal','MAST_Normal','EPI_Normal',
        'ENDO_Normal','DT_Normal')
sigScores<-sigScores[,c(adpkd,control)]

#创建X轴细胞所属样本类型标签
annotation_col<-data.frame(Group=factor(c(rep('ADPKD',9),rep('Normal',9))),
                           Type=factor(rep(c('PT','Tcell','MNP',
        'FIB','Bcell','MAST','EPI',
        'ENDO','DT'),2)))
rownames(annotation_col)<-colnames(sigScores)

#标签颜色
library(ggsci)
ann_colors<-list(Group=c('ADPKD'='#f47720','Normal'='#459943'),
                       Type=c('PT'=pal_npg()(9)[1],'Tcell'=pal_npg()(9)[2],
                             'MNP'=pal_npg()(9)[3],'FIB'=pal_npg()(9)[4],
                             'Bcell'=pal_npg()(9)[5],'MAST'=pal_npg()(9)[6],
                             'EPI'=pal_npg()(9)[7],'ENDO'=pal_npg()(9)[8],
                             'DT'=pal_npg()(9)[9]))

#   Fig: VISION结果热图
#   Fig: Hallmark
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

#   Fig: Hallmark
png(filename = "ADPKD_VISION_hallmark_all.png",width = 230, height = 280, 
    units = "mm", res = 300)
p
dev.off()
