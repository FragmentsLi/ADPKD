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
#   4. Pathway enrichment
####################################################################
library(ggplot2)
library(ggsci)
library(VISION)
Tcell <- subset(Tcell, subset = batch == "ADPKD")
Tcell <- subset(Tcell, subset = subtype != "NK")


#hallmark
vision.obj <- Vision(Tcell,
              signatures = c(".../msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

#GOBP
vision.obj <- Vision(Tcell,
              signatures = c(".../msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
group.by <- "subtype" 
#mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
sigScores <- expm1(x = sigScores) 
sigScores <- aggregate(sigScores, by=list(Tcell@meta.data[[group.by]]), FUN="mean") 
rownames(sigScores) <- sigScores$Group.1 
sigScores <- t(sigScores[,-1]) 

####################################################################
#   5. Communication analysis
####################################################################
library(CellChat)
Tcell=subset(x = Tcell, subset = batch == "ADPKD")
data.input <- GetAssayData(Tcell, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Tcell)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#CellChat
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

#Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,population.size=T,type = "truncatedMean",trim = 0.001)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

####################################################################
#   6. Pseudotime trajectory analysis
####################################################################
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggsci)

Tcell <- subset(Tcell,idents = c('CD4+TNF+ C2','CD4+RTKN2+ C3'))
Tcell <- subset(x = Tcell, subset = batch == "ADPKD")

Tcell <- NormalizeData(Tcell)
set.seed(1218)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, features = row.names(Tcell),verbose = FALSE)
set.seed(1218)
Tcell <- RunPCA(Tcell, pc.genes = VariableFeatures(object = Tcell))
Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:20)
Tcell <- FindNeighbors(Tcell, reduction = "harmony", dims = 1:20)
Tcell <- RunTSNE(Tcell,dims=1:20)

data <- GetAssayData(Tcell, assay = 'RNA', slot = 'counts')
cell_metadata <- Tcell@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, preprocess_method = "PCA")


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Tcell, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

psofa<-as.numeric(labdata$PSOFA)
psofa1<-unlist(lapply(psofa,function(x){
  if(x<2) return('PSOFA<2')
  else return('PSOFA>=2')
}))

#assign paritions
re.partition<-c(rep(1,length(cds@colData@rownames)))
names(re.partition)<-cds@colData@rownames
re.partition<-as.factor(re.partition)
cds@clusters$UMAP$partitions<-re.partition

#assign the cluster info
cds@clusters$UMAP$clusters<-as.character(Tcell@meta.data$subtype)

cds <- learn_graph(cds,use_partition=FALSE)
