library(Biobase)
library(Seurat) #Seurat4
library(Libra)
# Libra for Seurat4 can be found at https://github.com/jumphone/public/tree/main/download

SEED=123 # We have tested other seed and got similar results.
pbmc=readRDS('../Cerebellum/fetel_cerebellum_SoupX_nf2_dc_ep.rds')

# because development stages (different weeks) are not replicates, we ramdomly assign replicate labels
set.seed(SEED)
pbmc$replicate=sample(c('S1','S2','S3','S4','S5'),ncol(pbmc),replace=T)

# because the function named "run_de" use "label" to guide the  DEG analysis, we assign cell types ('TCP' v.s. "nonTCP") to the variable named "label" but not "cell_type".
pbmc$cell_type=rep('All',ncol(pbmc))
pbmc$label=rep('nonTCP',length(pbmc$cell_type))
pbmc$label[which(pbmc$celltype=='02xTCP')]='TCP'
pbmc$label=factor(pbmc$label,levels = c("TCP","nonTCP"))

library(Libra)
library(dplyr)
DE_LST=list()

DE_LST$id1 = run_de(pbmc,n_threads = 20,de_family = "pseudobulk",de_method = "edgeR", de_type='LRT')
DE_LST$id2 = run_de(pbmc,n_threads = 20,de_family = "pseudobulk",de_method = "edgeR", de_type='QLF')
DE_LST$id3 = run_de(pbmc,n_threads = 20,de_family = "pseudobulk",de_method = "DESeq2", de_type='LRT')
DE_LST$id4 = run_de(pbmc,n_threads = 20,de_family = "pseudobulk",de_method = "DESeq2", de_type='Wald')
DE_LST$id5 = run_de(pbmc,n_threads = 20,de_family = "pseudobulk",de_method = "limma", de_type='trend')
DE_LST$id6 = run_de(pbmc,n_threads = 20,de_family = "pseudobulk",de_method = "limma", de_type='voom')

###################################
saveRDS(DE_LST,file='DE_LST_ALL.rds')

OUT=list()
OUT$edger_lrt=DE_LST$id1
OUT$edger_qlf=DE_LST$id2
OUT$deseq2_lrt=DE_LST$id3
OUT$deseq2_wald=DE_LST$id4
OUT$limma_trend=DE_LST$id5
OUT$limma_voom=DE_LST$id6

saveRDS(OUT,file='DE_OUT_ALL.rds')
###################################

LGFD_CTNNB1=c()
PADJ_CTNNB1=c()
LGFD_SOX11=c()
PADJ_SOX11=c()
LGFD_HNRNPH1=c()
PADJ_HNRNPH1=c()

i=1
NAME=names(OUT)
while(i<=length(OUT)){
  this_out=OUT[[i]]
  this_name=NAME[i]
  this_lgfc=this_out$avg_logFC
  this_padj=this_out$p_val_adj
  this_ctnnb1=which(this_out$gene=='CTNNB1')
  this_sox11=which(this_out$gene=='SOX11')
  this_hnrnph1=which(this_out$gene=='HNRNPH1')
  LGFD_CTNNB1=c(LGFD_CTNNB1, this_lgfc[this_ctnnb1])
  PADJ_CTNNB1=c(PADJ_CTNNB1, this_padj[this_ctnnb1])
  LGFD_SOX11=c(LGFD_SOX11, this_lgfc[this_sox11])
  PADJ_SOX11=c(PADJ_SOX11, this_padj[this_sox11])
  LGFD_HNRNPH1=c(LGFD_HNRNPH1, this_lgfc[this_hnrnph1])
  PADJ_HNRNPH1=c(PADJ_HNRNPH1, this_padj[this_hnrnph1])
  i=i+1}


LGFD_MAT=rbind(LGFD_CTNNB1, LGFD_SOX11, LGFD_HNRNPH1)
PADJ_MAT=rbind(PADJ_CTNNB1, PADJ_SOX11, PADJ_HNRNPH1)
colnames(LGFD_MAT)=NAME
colnames(PADJ_MAT)=NAME
rownames(LGFD_MAT)=c('CTNNB1','SOX11','HNRNPH1')
rownames(PADJ_MAT)=c('CTNNB1','SOX11','HNRNPH1')

###########################################

library('ComplexHeatmap')
library('circlize')

col_fun =colorRamp2(c(-2, -1, 0, 1, 2 ), c('royalblue3','royalblue1','grey90','indianred1','red1'))
o.mat=as.matrix(LGFD_MAT)
ht1=Heatmap(o.mat,row_title='',name="lgfc",cluster_rows=F,cluster_columns=F,
            show_column_dend = F, show_row_dend = F,
            show_column_names=T, show_row_names=T,
            col=col_fun, border = TRUE,
            row_names_side='right',
            row_names_gp = gpar(fontsize = 10, lineheight=NULL),
            column_names_gp = gpar(fontsize = 10, lineheight=NULL),
            row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
            column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }
)

print(ht1)


col_fun =colorRamp2(c(0.001, 0.01, 0.05, 0.1, 1 ), c('red1','indianred1','grey90','grey80','royalblue1'))
o.mat=as.matrix(PADJ_MAT)
ht2=Heatmap(o.mat,row_title='',name="padj",cluster_rows=F,cluster_columns=F,
            show_column_dend = F, show_row_dend = F,
            show_column_names=T, show_row_names=T,
            col=col_fun, border = TRUE,
            row_names_side='right',
            row_names_gp = gpar(fontsize = 10, lineheight=NULL),
            column_names_gp = gpar(fontsize = 10, lineheight=NULL),
            row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
            column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.3f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }
)

print(ht2)


pdf('LGFD_PADJ.pdf',width=5,height=3)
print(ht1)
print(ht2)
dev.off()



write.table(OUT$edger_lrt,'OUT_edger_lrt.txt',sep='\t',row.names=T,col.names=T,quote=F)