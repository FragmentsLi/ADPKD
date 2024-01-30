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

