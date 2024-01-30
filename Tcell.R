####################################################################
#   1. Importing Data
####################################################################
library(Seurat)
library(ggplot2)
load('adpkd_all.RData')

Tcell <- subset(adpkd_qc,idents = c("Tcell"))
colnames(Tcell@meta.data)[5]<-'CellType.0.7'

####################################################################
#   2. Preprocessing
####################################################################
Tcell <- NormalizeData(Tcell)
set.seed(1218)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, features = row.names(Tcell),verbose = FALSE)
set.seed(1218)
Tcell <- RunPCA(Tcell, pc.genes = VariableFeatures(object = Tcell))

#  harmony
library(harmony)
set.seed(1218)
Tcell <- RunHarmony(Tcell,"batch", plot_convergence = TRUE)
Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:20)
Tcell <- FindNeighbors(Tcell, reduction = "harmony", dims = 1:20)

library(clustree)
set.seed(1218)
Tcell <- FindClusters(
  object = Tcell,
  resolution = c(seq(.1,1.6,.2))
)

####################################################################
#   3. Annotation
####################################################################
Idents(Tcell) <- 'RNA_snn_res.0.1'
#remove cluster 7
Tcell<-subset(Tcell,idents=c(0:6,8))
type<-c('CD4+IL7R+ C1','CD8+NKG7+ C1','CD8+CCL4L2+ C2','CD4+TNF+ C2','CD8+MGP+ C3',#0,1,2,3,4
        'NK','CD4+RTKN2+ C3','CD8+STMN1+ C4')#5,6,8
names(type) <- levels(Tcell) 
Tcell <- RenameIdents(Tcell, type)
Idents(Tcell)<-ordered(Idents(Tcell),
                       levels=c('CD4+IL7R+ C1','CD4+TNF+ C2','CD4+RTKN2+ C3',
                                'CD8+NKG7+ C1','CD8+CCL4L2+ C2','CD8+MGP+ C3','CD8+STMN1+ C4','NK'))
Tcell$subtype <- Idents(Tcell)
####################################################################
#   4. Annotation
####################################################################

####################################################################
#   5. Communication analysis
####################################################################
