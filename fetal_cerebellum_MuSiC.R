

library(MuSiC)
library(Biobase)
library(Seurat)


pbmc=readRDS('../Cerebellum/pbmc_sub.rds')##Only use neuronal cells for MuSiC

TMP=read.table('MB1_RNA_for_deconv_all.csv',sep=',',row.names=NULL,header=T)
TMP=TMP[which(TMP[,1] %in% names(which(table(TMP[,1])==1))),]
BULK=TMP[,c(2:ncol(TMP))]
rownames(BULK)=TMP[,1]
ANNO=t(read.table('MB1_RNA_for_deconv_all.csv.type',sep=','))[,1][2:(ncol(BULK)+1)]
colnames(BULK)=paste0(ANNO,'_',colnames(BULK))




set.seed(123)
USED_INDEX=sample(c(1:ncol(pbmc)),10000)
MAT=pbmc[['RNA']]@data[,USED_INDEX]
#MAT=MAT[which(rownames(MAT) %in% rownames(BULK)),]
MAT=as.matrix(MAT)
PD= cbind(as.character(pbmc$type)[USED_INDEX], as.character(pbmc$orig.ident)[USED_INDEX])
rownames(PD)=colnames(MAT)
colnames(PD)=c('type','sample')
PD=as.data.frame(PD)
#PD <- new("AnnotatedDataFrame", data=PD)
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=MAT),colData=DataFrame(type=PD$type, sample=PD$sample),,metadata=PD)
ES.REF=sce


MAT=as.matrix(BULK)
GS= ExpressionSet(assayData=MAT)
PD= cbind(colnames(MAT),colnames(MAT))
rownames(PD)=colnames(MAT)
colnames(PD)=c('type1','type2')
PD=as.data.frame(PD)
PD <- new("AnnotatedDataFrame", data=PD)
ES <- ExpressionSet(assayData=MAT,phenoData=PD)
ES.BULK=ES


set.seed(1234)
Est.bulk = music_prop(bulk.mtx = exprs(ES.BULK), sc.sce = ES.REF, clusters = 'type',samples='sample')

saveRDS(Est.bulk,'MUSIC_result.rds')

MUSIC=Est.bulk$Est.prop.weighted
MUSIC=MUSIC[,order(colnames(MUSIC))]

write.table(MUSIC,file='MUSIC_prop.txt',sep='\t',quote=F,row.names=T,col.names=T)








