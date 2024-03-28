library(dplyr)
library(Seurat)
library(patchwork)
library(SoupX)

#remove the ambient RNA form the cells using SoupX for each sample
sc= load10X('E:/..../FB20191')
sc = autoEstCont(sc)
out = adjustCounts(sc)
pbmc <- CreateSeuratObject(counts = out, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(pbmc), 20)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")

## generate or use the nf2 file for each sample
library(DropletQC)

nf2 <- nuclear_fraction_annotation(
  annotation_path = "/..../hg38_genes.gtf",
  bam ="/..../possorted_genome_bam.bam",
  barcodes = "/..../barcodes.tsv.gz",
  tiles = 1, cores = 1, verbose = FALSE)
saveRDS(nf2, file = "/..../FB20191_nf2.rds")

nf2=readRDS("E:/..../FB20191_nf2.rds")
pbmc<-AddMetaData(pbmc, nf2)

input=pbmc@meta.data[,c("nuclear_fraction","nCount_RNA")]
input2 = identify_empty_drops(nf_umi = input)
input2$seurat_clusters <- pbmc$seurat_clusters
input2.dc <- identify_damaged_cells(input2, verbose=FALSE)
pbmc<-AddMetaData(pbmc, input2.dc[[1]])
head(pbmc)
DimPlot(pbmc, group.by = "cell_status")

saveRDS(pbmc, file = "E:/..../FB20191_SoupX_nf2_dc.rds")


#Merge these samples together

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

PCW8 = readRDS("/..../FB15_SoupX_nf2_dc.rds")
PCW9 = readRDS("/..../FB13_SoupX_nf2_dc.rds")
PCW12 = readRDS("/..../FB20195_SoupX_nf2_dc.rds")
PCW13 = readRDS("/..../FB20197_SoupX_nf2_dc.rds")
PCW14 = readRDS("/...../FB20192_SoupX_nf2_dc.rds")
PCW15 = readRDS("/..../FB20191_SoupX_nf2_dc.rds")
PCW16 = readRDS("/..../FB20193_SoupX_nf2_dc.rds")
PCW17 = readRDS("/..../FB20194_SoupX_nf2_dc.rds")

colnames(PCW8)<-paste0("PCW8_",colnames(PCW8))
colnames(PCW9)<-paste0("PCW9_",colnames(PCW9))
colnames(PCW12)<-paste0("PCW12_",colnames(PCW12))
colnames(PCW13)<-paste0("PCW13_",colnames(PCW13))
colnames(PCW14)<-paste0("PCW14_",colnames(PCW14))
colnames(PCW15)<-paste0("PCW15_",colnames(PCW15))
colnames(PCW16)<-paste0("PCW16_",colnames(PCW16))
colnames(PCW17)<-paste0("PCW17_",colnames(PCW17))


pbmc = merge(x = PCW8, y = c(PCW9, PCW12, PCW13, PCW14, PCW15, PCW16, PCW17))

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ElbowPlot(pbmc)

BATCH=c(rep('PCW8',ncol(PCW8)),
        rep('PCW9',ncol(PCW9)),
        rep('PCW12',ncol(PCW12)),
        rep('PCW13',ncol(PCW13)),
        rep('PCW14',ncol(PCW14)),
        rep('PCW15',ncol(PCW15)),
        rep('PCW16',ncol(PCW16)), 
        rep('PCW17',ncol(PCW17)))
pbmc$batch=BATCH

pbmc <- FindNeighbors(pbmc, dims = 1:23)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:23)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", raster=FALSE)
new.cluster.ids <- c("GCP-1", "UBC", "GCP-2", "GCP-pro", "Granule_cell-1", "Diff.UBC", "Cycling_GCP-pro", "TCP", "NSC", "Microglia", "Dev.Purkinje", "Meninges", "Microglia", "UBC-pro", "GABA_interneuron", "Granule_cell-2", "NSC", "Brainstem", "Purkinje", "Meninges", "Brainstem", "Endothelial", "Unkown", "UBC", "Monocyte", "Glial.progenitor", "Pericyte", "NSC", "Dev.Purkinje", "Vascular.cell", "Red_blood_cell", "Ependymal", "Oligodendrocytes")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype=pbmc@active.ident
DimPlot(pbmc, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.5)
DimPlot(pbmc, group.by = "cell_status", raster=FALSE)
table(pbmc$celltype, pbmc$cell_status)
saveRDS(pbmc, file = "/..../fetel_cerebellum_SoupX_nf2_dc_ep.rds")