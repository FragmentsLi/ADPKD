####################################################################
#   1. Microarray dataset preprocessing
####################################################################
library(affy)
library(stringr)

setwd('.../GSE7869')
data.raw <- ReadAffy()
sampleNames(data.raw)

# RMA
data.rma <- rma(data.raw)
data<-exprs(data.rma)
colnames(data)<-unlist(lapply(colnames(data),function(x){return(str_sub(x,1,9))}))
#21 samples, 54,675 probes

map<-read.delim('GPL570-55999.txt',header=T,skip=16)
map<-map[,c('ID','Gene.Symbol')]

# Probe to Gene symbol, get max median probe corresponding to one gene symbol
map2<-data.frame(matrix(nrow=0,ncol=2))
colnames(map2)<-c('id','genesymbol')

for(i in 1:nrow(map)){
  if(str_detect(map[i,2],'///')){
    genes<-str_split(map[i,2],' /// ')[[1]]
    for(j in 1:length(genes)){
      map2<-rbind(map2,data.frame(id=map[i,1],genesymbol=genes[j]))}}
  else
    map2<-rbind(map2,data.frame(id=map[i,1],genesymbol=map[i,2]))}

data2<-data.frame(matrix(nrow=0,ncol=ncol(data)))
colnames(data2)<-colnames(data)

for(gene in unique(map2$genesymbol)){
  tmp<-subset(map2,map2$genesymbol==gene)
  if(gene=='')
    next
  if(nrow(tmp)==1){
    data2<-rbind(data2,data[tmp[1,1],])
    rownames(data2)[nrow(data2)]=gene
  }
  else{
    median_list<-rowMedians(as.matrix(data[intersect(tmp$id,rownames(data)),]))
    data2<-rbind(data2,data[names(median_list[order(median_list,decreasing=T)])[1],])
    rownames(data2)[nrow(data2)]=gene}}

write.csv(data2,'GSE7869_rma.csv',row.names=T,quote=F)

####################################################################
#   2. CIBERSORTx
####################################################################

library(Seurat)
load('adpkd_all.RData')
source('.../DMC-master/R/create_sig_matrix.R')

# Get normalized matrix
df <- as.data.frame(GetAssayData(adpkd_qc, assay = "RNA", slot = "counts"))
# Retain HVGs
df<-df[VariableFeatures(object = adpkd_qc),]

pheno<-data.frame(cell_type=Idents(adpkd_qc))

adpkd_qc=0

# Create signature matrix
options(mc.cores=10)
sig_matrix<-create_sig_matrix(df,pheno)

# Microarray data
library(IOBR)
mixture_file<-read.csv('GSE7869_rma.csv',row.names=1,check.names=F)

options(mc.cores=10)

cibersortx_adpkd<-CIBERSORT(
  sig_matrix=sig_matrix,
  mixture_file=mixture_file,
  perm=1000,absolute=FALSE)
