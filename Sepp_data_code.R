
##Seurat v4
library(dplyr)
library(Seurat)
library(patchwork)
pbmc = readRDS("/Volumes/Zaili/cellranger3.1_gh38/Sepp.rds")

DATA=pbmc[['RNA']]@data

library(clusterProfiler)
library(org.Hs.eg.db)
ENSG=rownames(DATA)
SYMBOL=mapIds(org.Hs.eg.db, ENSG, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

USED=which(!is.na(SYMBOL))

DATA=DATA[USED,]
rownames(DATA)=SYMBOL[USED]

###DATA=pbmc@assays$RNA@data

RDATA=DATA[,sample(1:ncol(DATA),120000)]
duplicated_rows <- RDATA[duplicated(rownames(RDATA)), ]
if (nrow(duplicated_rows) > 0) {
  rownames(RDATA) <- make.names(rownames(RDATA), unique=TRUE)
}

pbmc <- CreateSeuratObject(counts = RDATA, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
#QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)
rm(DATA)
rm(RDATA)
gc()

###pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:23)
pbmc <- FindClusters(pbmc, resolution = 0.8)

pbmc <- RunUMAP(pbmc, dims = 1:14)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = T)
FeaturePlot(pbmc, features = c("NES", "HES5", "SOX11", "HNRNPH1"), cols = c('azure3',"darkgoldenrod","red", "red3"))


new.cluster.ids <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19","c20", "c21", "c22", "c23", "c24", "c25", "c26", "c27", "c28", "c29","c30", "c31", "c32", "c33", "c34", "c35","c36", "c37", "c38", "c39", "c40", "c41")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

pbmc@meta.data$celltype=pbmc@active.ident
saveRDS(pbmc, file = "/Volumes/Zaili/cellranger3.1_gh38/Sepp_120k_cells.rds")




##Seurat 5
##rliger 1.0.0
library(rliger)
library(Seurat)
library(SeuratWrappers)
d2 = readRDS("D:/cellranger3.1_gh38/Sepp_120k_cells.rds")
DimPlot(d2, reduction = "umap", label = TRUE, pt.size = 0.5)
d2@meta.data$celltype=d2@active.ident
d1 = readRDS("D:/zaili_old_driver/Data_all/Data/single_cells/TOTAL/Refere/LIGER/Human23_1final.rds")

d3 <- merge(d1,d2) ##  62415 122146
g1 <- rownames(d1@assays$RNA@data);g2 <- rownames(d2@assays$RNA@data)
g12 <- intersect(g1,g2)
#d3 <- merge(d1,d2,add.cell.ids=c('old','new'))
d3@meta.data$original_dataset <- ifelse(d3@meta.data$orig.ident %in% d1@meta.data$orig.ident,'scRNA','snRNA')
w1 <- which(d3@meta.data$original_dataset=='old')
w2 <- which(d3@meta.data$original_dataset=='new')
d3@meta.data$CellType <- d3@meta.data$celltype
d3@meta.data$CellType[w2] <- d3@meta.data$fig_cell_type[w2]
##
pbmcsca <- d3
rm(d1)
rm(d2)
rm(d3)
gc()

pbmcsca <- NormalizeData(pbmcsca)
pbmcsca <- FindVariableFeatures(pbmcsca)
pbmcsca <- ScaleData(pbmcsca, split.by = "original_dataset", do.center = FALSE)
pbmcsca <- RunOptimizeALS(pbmcsca, k = 16, lambda = 5, split.by = "original_dataset")
pbmcsca <- RunQuantileNorm(pbmcsca, split.by = "original_dataset")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs

pbmcsca <- FindNeighbors(pbmcsca, reduction = "iNMF", dims = 1:15)
pbmcsca <- FindClusters(pbmcsca, resolution = 0.4)
# Dimensional reduction and plotting
pbmcsca <- RunUMAP(pbmcsca, dims = 1:ncol(pbmcsca[["iNMF"]]), reduction = "iNMF")
DimPlot(pbmcsca, group.by = c("original_dataset", "CellType"), ncol = 2)
##
pbmcsca.liger <- seuratToLiger(pbmcsca,combined.seurat = T,
                               meta.var = 'original_dataset')
cluster1 <- factor(pbmcsca@meta.data$CellType[which(pbmcsca@meta.data$original_dataset=='scRNA')])
cluster2 <- factor(pbmcsca@meta.data$CellType[which(pbmcsca@meta.data$original_dataset=='snRNA')])
names(cluster1) <- rownames(pbmcsca@meta.data)[which(pbmcsca@meta.data$original_dataset=='scRNA')]
names(cluster2) <- rownames(pbmcsca@meta.data)[which(pbmcsca@meta.data$original_dataset=='snRNA')]
par(mar=c(1,8,1,8));par(xpd=T)
makeRiverplot(pbmcsca.liger,
              cluster1 = cluster1,
              cluster2 = cluster2,river.yscale = 7,label.cex = 1, min.frac = 0.1)

saveRDS(pbmcsca, file = "D:/cellranger3.1_gh38/Sepp/Sepp_pbmcsca.rds")



