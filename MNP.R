####################################################################
#   1. Importing Data
####################################################################
library(Seurat)
library(ggplot2)
load('adpkd_all.RData')

MNP <- subset(adpkd_qc,idents = c("MNP"))
colnames(MNP@meta.data)[5]<-'CellType.0.7'

####################################################################
#   2. Preprocessing
####################################################################
MNP <- NormalizeData(MNP)
set.seed(1218)
MNP <- FindVariableFeatures(MNP, selection.method = "vst", nfeatures = 2000)
MNP <- ScaleData(MNP, features = row.names(MNP),verbose = FALSE)
set.seed(1218)
MNP <- RunPCA(MNP, pc.genes = VariableFeatures(object = MNP))

#  harmony
library(harmony)
set.seed(1218)
MNP <- RunHarmony(MNP,"batch", plot_convergence = TRUE)
MNP <- RunUMAP(MNP, reduction = "harmony", dims = 1:20)
MNP <- FindNeighbors(MNP, reduction = "harmony", dims = 1:20)
MNP <- RunTSNE(MNP,dims=1:20)
library(clustree)
set.seed(1218)
MNP <- FindClusters(
  object = MNP,
  resolution = c(seq(.1,1.6,.2))
)

####################################################################
#   3. Annotation
####################################################################
Idents(MNP) <- 'RNA_snn_res.0.3'
type<-c('S100A8+ cMono','SELENOP+ M2_C1','IL1B+ M1_C1','GPX3+ M2_C2','SELENOP+ M2_C1',#0,1,2,3,4
        'CCL5+ M1_C2','APOC1+ M2_C3','APOC1+ M2_C3','GZMB+ pDC_C1','S100B+ pDC_C2',#5,6,7,8,9
        'CCL5+ M1_C2')#10

names(type) <- levels(MNP) 
MNP <- RenameIdents(MNP, type) 
Idents(MNP)<-ordered(Idents(MNP),levels=c('IL1B+ M1_C1','CCL5+ M1_C2','SELENOP+ M2_C1',
                                          'GPX3+ M2_C2','APOC1+ M2_C3','S100A8+ cMono',
                                          'GZMB+ pDC_C1','S100B+ pDC_C2'))

####################################################################
#   4. M1/M2 score
####################################################################
M1<-c('CCR7','IL2RA','IL15RA','IL7R',#Membrane receptors 
      'CXCL11','CCL19','CXCL10','CXCL9','TNF','CCL5','CCL15','IL12B','IL15','TRAIL','IL6','CCL20','PBEF1','ECGF1',#Cytokines and chemokines 
      'BCL2A1','FAS','BIRC3','GADD45G','HSXIAPAF1',#Apoptosis-related genes 
      'SLC7A5','SLC21A15','SLC2A6','SLC31A2',#Solute carriers
      'INDO','PLA1A','OASL','CHI3L2','HSD11B1','AK3','SPHK1','PFKFB3','PSME2','PFKP','PSMB9','PSMA2','OAS2',#Enzymes
      'PTX3','CSPG2','APOL3','IGFBP4','APOL1','PDGFA','EDN1','APOL2','INHBA','APOL6',#Extracelullar mediators
      'HESX1','IRF1','ATF3','IRF7'#DNA-binding factors 
     )

M2<-c('GPR86','P2RY5','TGFBR2','HRH1','TLR5','DCL-1','MSR1','CXCR4','DECTIN1','P2RY14','DCSIGN','CLECSF13','MS4A6A','CD36','MS4A4A','MRC1',##Membrane receptors
      'IGF1','CCL23','CCL18','CCL13',#Cytokines and chemokines
      'SLC21A9','SLC4A7','SLC38A6',#Solute carriers
      'CTSC','HEXB','LIPA','ADK','HNMT','TPST2','CERK','HS3ST2','LTA4H','CA2','ALOX15','HS3ST1',#Enzymes
      'TGFBI','SEPP1','CHN2','FN1','FGL2',#Extracelullar mediators
      'GAS7','EGR2','MAF'#DNA-binding factors 
     )

M1_mean<-colMeans(MNP[['RNA']]@scale.data[intersect(rownames(MNP[['RNA']]@scale.data),M1),])
M2_mean<-colMeans(MNP[['RNA']]@scale.data[intersect(rownames(MNP[['RNA']]@scale.data),M2),])

MNP@meta.data<-data.frame(MNP@meta.data,M1=M1_mean[rownames(MNP@meta.data)],M2=M2_mean[rownames(MNP@meta.data)])

####################################################################
#   5. Pathway enrichment
####################################################################
library(VISION)
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
                             
                             
