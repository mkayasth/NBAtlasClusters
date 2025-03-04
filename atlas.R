library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)

# needs to be set for large dataset analysis
options(future.globals.maxSize = 6 * 1024^3) # setting max to 6 GB.

# loading the rds file from NBAtlas -- lightweight version.
atlas <- readRDS("NBAtlas/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds")

# first, lets run supervised data cluster from the scVI_umap that came along with the rds file.
dim(Embeddings(atlas, reduction = "scvi_umap"))
atlas <- FindNeighbors(atlas, reduction = "scvi_umap", dims = 1:2)  # scVI UMAP has 2 dimensions for 50k cells.
atlas <- FindClusters(atlas, resolution = 0.5)

table(atlas$orig.ident)
table(atlas$Cell_type)
# orig.ident and cell_type data match exactly in the .rds file (original annotation and annotation from scVI umap are the same).

# dim plot for the scvi umap obtained from the paper itself.
plot_supervised <- DimPlot(atlas, reduction = "scvi_umap", group.by = "orig.ident")
plot_supervised

# Now, running PCA, umap and seeing if the clusters still match.
atlas <- readRDS("NBAtlas/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds")
atlas2 <- NormalizeData(atlas)
atlas2 <- FindVariableFeatures(atlas2)

# unsupervised clustering of the cells.
atlas2 <- ScaleData(atlas2)
atlas2 <- RunPCA(atlas2)

dim(Embeddings(atlas2, reduction = "pca"))
atlas2 <- FindNeighbors(atlas2, dims = 1:50)
atlas2 <- FindClusters(atlas2, resolution = 0.5)
atlas2 <- RunUMAP(atlas2, dims = 1:50, return.model = T)


# Visualize UMAP based on supervised clustering.
plot_unsupervised <- DimPlot(atlas2, reduction = "umap", group.by = "orig.ident")

# clusters as per umap (unsupervised).
DimPlot(atlas2, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()

# unsupervised clusters but data labeled later from orig.identities.
plot_unsupervised


plot_supervised + plot_unsupervised

# Maybe batch correction is the issue here. We have data from 7 different studies.
atlas <- readRDS("NBAtlas/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds")
atlas2 <- NormalizeData(atlas)
atlas2 <- FindVariableFeatures(atlas2)
atlas2 <- ScaleData(atlas2)
atlas2 <- RunPCA(atlas2)

atlas2 <- RunHarmony(atlas2, group.by.vars = "Study")

# number of dimensions in this reduction.
atlas2@reductions$harmony


atlas2 <- FindNeighbors(atlas2, reduction = "harmony", dims = 1:50)
atlas2 <- FindClusters(atlas2, resolution = 0.5)

atlas2 <- RunUMAP(atlas2, reduction = "harmony", dims = 1:50)


plot_unsupervised_batchEffect <- DimPlot(atlas2, reduction = "umap", group.by = "orig.ident")
plot_unsupervised_batchEffect

# self label based on cluster.
DimPlot(atlas2, reduction = "umap", group.by = "seurat_clusters")

plot_supervised + plot_unsupervised_batchEffect
