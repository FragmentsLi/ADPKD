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

# Creating Figure 4B
library(ggsci)
library(reshape2)
library(ggplot2)

colnames(cibersortx_adpkd)<-c('MNP','Tcell','MAST','Bcell','EPI','PT','ENDO','FIB','DT')
cibersortx_adpkd<-as.matrix(cibersortx_adpkd[,1:9])
cibersortx_adpkd<-melt(cibersortx_adpkd,id=rownames(cibersortx_adpkd))
colnames(cibersortx_adpkd)<-c('Sample','Group','Cell')

data.anno<-data.frame(Sample=unique(cibersortx_adpkd$Sample),
                      label=c(rep('Small Cyst',5),rep('Medium Cyst',5),rep('Large Cyst',3),
                             rep('Minimal Cyst',5),rep('Normal',3)))

df.plot<-plyr::join(cibersortx_adpkd,data.anno)

df.plot$Group<-ordered(df.plot$Group,
  levels=c('PT','Tcell','MNP','FIB',
           'Bcell','MAST','EPI','ENDO','DT'))
df.plot$Sample<-ordered(df.plot$Sample,
  levels=c("GSM190876","GSM190877","GSM190878","GSM190871","GSM190872","GSM190873",
           "GSM190874","GSM190875","GSM190858","GSM190859","GSM190860","GSM190861",
           "GSM190862","GSM190863","GSM190864","GSM190865","GSM190866","GSM190867",
           "GSM190868","GSM190869","GSM190870"))

df.plot$label<-ordered(df.plot$label,
  levels=c('Normal','Minimal Cyst','Small Cyst','Medium Cyst','Large Cyst'))

p<-ggplot(df.plot,aes(x=Sample,y=Cell,fill=Group))+
geom_bar(stat="identity", position="stack", width =0.8 ,aes(fill=Group))+
labs(x = "",y = "Proportion of Cells")+theme_bw()+
theme(axis.ticks.length=unit(0.1,'cm'))+
guides(fill=guide_legend(title=NULL))+
scale_fill_manual(values=pal_npg()(9))+ 
facet_grid(~label,scales="free_x",space="free")+theme_classic()+
theme(plot.title=element_text(size = 8,hjust=0.5,face='bold'),
      legend.text=element_text(size=8),
          legend.title=element_text(size=8),
     legend.position='bottom') + 
theme(axis.title = element_text(size = 8,face='bold'),
          axis.text.y=element_text(size=8,color = "black"),
          axis.text.x = element_text(size = 8,hjust=1,vjust=0.5,angle=90,color = "black"), 
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1)) 

png(filename = "ADPKD_CIBERSORTx_HVGs.png",width = 155, height = 85, 
    units = "mm", res = 300)
p
dev.off()
