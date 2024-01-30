library(CellChat)
library(Seurat)
load('adpkd_all.RData')

####################################################################
#   1. CellChat for ADPKD
####################################################################

# ADPKD sample
adpkd_qc<-subset(x=adpkd_qc,batch=='ADPKD')
data.input <- GetAssayData(adpkd_qc, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(adpkd_qc)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#CellChat
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# LR database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

#Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Calculation
cellchat <- computeCommunProb(cellchat,population.size=T,type = "truncatedMean",trim = 0.001)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

data.input=0
meta=0
cellchat_adpkd=cellchat

####################################################################
#   2. CellChat for Normal
####################################################################

# Normal sample
normal<-subset(x=adpkd_qc,batch=='Normal')
data.input <- GetAssayData(normal, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(normal)
#labels <- adpkd_qc@meta.data$subtype
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#CellChat
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# LR database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

#Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Calculation
cellchat <- computeCommunProb(cellchat,population.size=T,type = "truncatedMean",trim = 0.001)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

data.input=0
meta=0
cellchat_normal=cellchat

####################################################################
#   3. Comparison
####################################################################
object.list <- list(ADPKD = cellchat_adpkd, Normal = cellchat_normal)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat_adpkd=0
cellchat_normal=0

# Creating Figure 1E
png(filename = "ADPKD_CellChat_diff_heatmap.png",width = 100, height = 80, 
    units = "mm", res = 300)
netVisual_heatmap(cellchat,comparison = c(2, 1), 
                  title.name="Interaction weights/strength",color.use=pal_npg()(9),
                  measure = "weight",
                  font.size=8,
                 font.size.title = 8)

dev.off()
