
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

data <- Read10X(data.dir = "/Users/sravanisaadhu/Downloads/filtered_feature_bc_matrix/")
data <- CreateSeuratObject(counts = data, project = "pbmc10k", min.cells = 3, min.features = 200)
data


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(data@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- NormalizeData(data)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the 50 most highly variable genes
top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot1,geom_text_repel = TRUE,points=top10, xnudge = 0, ynudge = 0)
plot1 + plot2

all.genes <- rownames(data)
# data <- ScaleData(data, features = VariableFeatures(data)) skylar's note: this did not keep
# relevant genes in scale.data
data <- ScaleData(data, features = all.genes)

data <- ScaleData(data, vars.to.regress = "percent.mt")

data <- RunPCA(data, features = VariableFeatures(object = data))
# Examine and visualize PCA results a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca") + NoLegend()

DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(data)

data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
head(Idents(data), 5)

data <- RunUMAP(data, dims = 1:10)
DimPlot(data, reduction = "umap")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(data, ident.1 = 2)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(data, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
data.markers <- FindAllMarkers(data, only.pos = TRUE)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(data, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(data, features = c("CAMK4", "ANK3"))
# you can plot raw counts as well
VlnPlot(data, features = c("CSGALNACT1", "PHACTR2"), slot = "counts", log = TRUE)

FeaturePlot(data, features = c("CAMK4", "FHIT", "LEF1", "PRKCA", "MAML2", "ANK3", "PLCL1", "IGF1R", "ATP10A"))

data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(data, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
DimPlot(data, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
