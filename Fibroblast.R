####################################################################
#   1. Importing Data
####################################################################
library(Seurat)
library(ggplot2)
library(harmony)
load('adpkd_all.RData')
fibroblast <- subset(adpkd_qc,idents = c("FIB"))
colnames(fibroblast@meta.data)[5]<-'CellType.0.7'

####################################################################
#   2. Preprocessing
####################################################################
fibroblast <- NormalizeData(fibroblast)
set.seed(1218)
fibroblast <- FindVariableFeatures(fibroblast, selection.method = "vst", nfeatures = 2000)
fibroblast <- ScaleData(fibroblast, features = row.names(fibroblast),verbose = FALSE)
set.seed(1218)
fibroblast <- RunPCA(fibroblast, pc.genes = VariableFeatures(object = fibroblast))

library(harmony)
set.seed(1218)
fibroblast <- RunHarmony(fibroblast,"batch", plot_convergence = TRUE)
fibroblast <- RunUMAP(fibroblast, reduction = "harmony", dims = 1:20)
fibroblast <- FindNeighbors(fibroblast, reduction = "harmony", dims = 1:20)
fibroblast <- RunTSNE(fibroblast,dims=1:20)

library(clustree)
set.seed(1218)
fibroblast <- FindClusters(
  object = fibroblast,
  resolution = c(seq(.1,1.6,.2))
)

####################################################################
#   3. Annotation
####################################################################
Idents(fibroblast) <- 'RNA_snn_res.0.1'
type<-c('ACTA2+ Myo','SYT1+ FIB1','IGF1+ FIB2','CCL4+ FIB3','GRM8+ FIB4',#0,1,2,3,4
        'SERPINE2+ FIB5','CFD+ FIB6','MME+ FIB7')#5,6,7

names(type) <- levels(fibroblast) 
fibroblast <- RenameIdents(fibroblast, type)
fibroblast$subtype <- Idents(fibroblast)

####################################################################
#   4. Pathway enrichment
####################################################################
library(VISION)
fibroblast <- subset(fibroblast, subset = batch == "ADPKD")

#hallmark
vision.obj <- Vision(fibroblast,
              signatures = c(".../msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

#GOBP
vision.obj <- Vision(fibroblast,
              signatures = c(".../msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt"),
              pool=F)
options(mc.cores=5)
vision.obj <- analyze(vision.obj)

####################################################################
#   5. Communication analysis
####################################################################
fibroblast<-subset(x=fibroblast,batch=='ADPKD')

data.input <- GetAssayData(fibroblast, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(fibroblast)
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

fibroblast <- subset(x = fibroblast, subset = batch == "ADPKD")
fibroblast

fibroblast <- NormalizeData(fibroblast)
set.seed(1218)
fibroblast <- FindVariableFeatures(fibroblast, selection.method = "vst", nfeatures = 2000)
fibroblast <- ScaleData(fibroblast, features = row.names(fibroblast),verbose = FALSE)
set.seed(1218)
fibroblast <- RunPCA(fibroblast, pc.genes = VariableFeatures(object = fibroblast))
fibroblast <- RunUMAP(fibroblast, reduction = "harmony", dims = 1:20)
fibroblast <- FindNeighbors(fibroblast, reduction = "harmony", dims = 1:20)
fibroblast <- RunTSNE(fibroblast,dims=1:20)

data <- GetAssayData(fibroblast, assay = 'RNA', slot = 'counts')
cell_metadata <- fibroblast@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(fibroblast, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

head(colData(cds))

#assign paritions
re.partition<-c(rep(1,length(cds@colData@rownames)))
names(re.partition)<-cds@colData@rownames
re.partition<-as.factor(re.partition)
cds@clusters$UMAP$partitions<-re.partition

#assign the cluster info
cds@clusters$UMAP$clusters<-as.character(fibroblast@meta.data$subtype)

#assign UMAP embeddings coordinate
#cds@int_colData@listData$reducedDims$UMAP<-fibroblast@reductions$umap@cell.embeddings

#monocle3
#cds <-cluster_cells(cds)
cds <- learn_graph(cds,use_partition=FALSE)

#ordercell
cds <- order_cells(cds,reduction_method='UMAP',
                   root_cells=rownames(subset(fibroblast@meta.data,fibroblast@meta.data$subtype=='IGF1+ FIB2')))

