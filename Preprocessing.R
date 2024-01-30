####################################################################
#   1. Importing Data
####################################################################
library(Seurat)

# ADPKD Sample ID
adpkd_id<-read.delim('SRR_Acc_List.txt',header=F)
adpkd_id<-adpkd_id[,1]

# Normal Sample ID
control_id<-c('GSM4145204','GSM4145205','GSM4145206')

# Path1
folders<-unlist(lapply(adpkd_id,function(x){
  paste0('/.../PRJNA679848/cellranger/',x,'/outs/filtered_feature_bc_matrix')
}))

# Path2
folders<-c(folders,unlist(lapply(control_id,function(x){
  paste0('/.../PRJNA679848/GSE131685_RAW/',x)
})))


# Cell Filtering: > 3 cells
adpkd_list<-lapply(c(1:length(folders)),function(x){
  CreateSeuratObject(count=Read10X(folders[x]),
                     project=c(adpkd_id,control_id)[x],min.cells=3)})

# Merge
adpkd_all<-merge(adpkd_list[[1]],
              y=c(adpkd_list[[2]],adpkd_list[[3]],adpkd_list[[4]],adpkd_list[[5]],adpkd_list[[6]],
                  adpkd_list[[7]],adpkd_list[[8]],adpkd_list[[9]],adpkd_list[[10]]),
              add.cell.ids=c(adpkd_id,control_id),
              project='ADPKD')
adpkd_list<-list()

# Forming 31,434 genes，111,601 cells

####################################################################
#   2. MT/RP/HB Thresholds
####################################################################
# MT
adpkd_all[['percent.mt']]<-PercentageFeatureSet(adpkd_all,pattern='^MT-')
# RP
adpkd_all[['percent.rp']]<-PercentageFeatureSet(adpkd_all,pattern='^RP')
# HB
HB.genes<-c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
adpkd_all[['percent.HB']]<-PercentageFeatureSet(adpkd_all,pattern=HB.genes)

#  Creating Sample ID
adpkd_all@meta.data$orig.ident<-as.factor(adpkd_all@meta.data$orig.ident)
adpkd_all@meta.data$orig.ident<-ordered(adpkd_all@meta.data$orig.ident,
                          levels=c("SRR13095044","SRR13095045","SRR13095046","SRR13095047","SRR13095048","SRR13095049","SRR13095050",
                                   "GSM4145204","GSM4145205","GSM4145206"))

adpkd_all@meta.data$orig.ident2<-unlist(lapply(as.character(adpkd_all@meta.data$orig.ident),function(x){
  if(strsplit(x,1,3)[[1]][1]=='SRR')
  return(paste0('A_S',which(names(table(adpkd_all@meta.data$orig.ident))==x)))
  else
  return(paste0('N_S',(which(names(table(adpkd_all@meta.data$orig.ident))==x)-7)))}))

adpkd_all@meta.data$orig.ident2<-ordered(adpkd_all@meta.data$orig.ident2,
                          levels=c("A_S1","A_S2","A_S3","A_S4","A_S5",
                                   "A_S6","A_S7","N_S1","N_S2","N_S3"))

Idents(adpkd_all)<-'orig.ident2'

####################################################################
#   3. Quality Control
####################################################################
adpkd_qc<-subset(adpkd_all,subset=nFeature_RNA>500
                 &nCount_RNA >500 &nCount_RNA <25000
                 &percent.mt<20
                 &percent.HB<5)
# Forming 31,434 genes，93,483 cells
adpkd_all<-c()

####################################################################
#   4. Preprocessing
####################################################################
# 1) Normalization
adpkd_qc<-NormalizeData(adpkd_qc)

# 2) Finding HVGs
adpkd_qc<-FindVariableFeatures(adpkd_qc,selection.method = "vst", nfeatures = 2000)

# 3) Scaling ALL
adpkd_qc<-ScaleData(adpkd_qc,features = row.names(adpkd_qc),verbose = FALSE)

# 4) PCA/UMAP
adpkd_qc<-RunPCA(adpkd_qc,pc.genes = VariableFeatures(object = adpkd_qc))
adpkd_qc<-FindNeighbors(adpkd_qc,dims=1:20)
adpkd_qc<-RunUMAP(adpkd_qc,dims=1:20)

# 5) Harmony
library(harmony)
library(ggplot2)
adpkd_id<-read.delim('SRR_Acc_List.txt',header=F)
adpkd_id<-adpkd_id[,1]
adpkd_qc@meta.data$batch<-unlist(lapply(adpkd_qc@meta.data$orig.ident,function(x){
  if(x %in% adpkd_id)
    return('ADPKD')
  else
    return('Normal')}))

# Batch correction
set.seed(1218)
adpkd_qc<-RunHarmony(adpkd_qc,"batch", plot_convergence = TRUE)
#UMAP
adpkd_qc<-FindNeighbors(adpkd_qc,reduction='harmony',dims=1:20)
adpkd_qc<-RunUMAP(adpkd_qc,reduction='harmony',dims=1:20)

# 6) Resolution
library(clustree)
adpkd_qc <- FindClusters(
  object = adpkd_qc,
  resolution = c(seq(.1,1.6,.2))
)
####################################################################
#   5. Saving
####################################################################
adpkd_qc@meta.data<-adpkd_qc@meta.data[,c(1:3,8,12,7)]
save(adpkd_qc,file='adpkd_all.RData')
