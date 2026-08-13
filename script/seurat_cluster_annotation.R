#!/usr/bin/env Rscript

suppressPackageStartupMessages({

library(Seurat)
library(SeuratObject)
library(zellkonverter)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(ggplot2)

})


args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
    stop("
Usage:

Rscript seurat_cluster_annotation.R input.h5ad output_dir

Example:

Rscript seurat_cluster_annotation.R sample.h5ad result/

")
}


input_h5ad <- args[1]
outdir <- args[2]


dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)


cat("Reading h5ad file...\n")


##################################################
# 1. Read h5ad
##################################################

sce <- readH5AD(
    input_h5ad
)


counts <- assay(
    sce,
    "X"
)


# gene x cell

seu <- CreateSeuratObject(
    counts = counts,
    project = "h5ad",
    min.cells = 0,
    min.features = 0
)


cat(
"Cells:",
ncol(seu),
"\nGenes:",
nrow(seu),
"\n"
)


##################################################
# 2. Normalization
##################################################

seu <- NormalizeData(
    seu,
    normalization.method="LogNormalize",
    scale.factor=10000
)


##################################################
# 3. Find variable genes
##################################################

seu <- FindVariableFeatures(
    seu,
    selection.method="vst",
    nfeatures=3000
)


##################################################
# 4. Scaling
##################################################

seu <- ScaleData(
    seu,
    features=VariableFeatures(seu)
)


##################################################
# 5. PCA
##################################################

seu <- RunPCA(
    seu,
    features=VariableFeatures(seu),
    npcs=50
)


##################################################
# 6. UMAP
##################################################

seu <- FindNeighbors(
    seu,
    dims=1:30
)


seu <- FindClusters(
    seu,
    resolution=0.5,
    algorithm=1
)


seu <- RunUMAP(
    seu,
    dims=1:30
)


##################################################
# 7. Marker genes
##################################################

markers <- FindAllMarkers(
    seu,
    only.pos=TRUE,
    min.pct=0.25,
    logfc.threshold=0.25
)


write.csv(
    markers,
    file=file.path(
        outdir,
        "cluster_markers.csv"
    ),
    row.names=FALSE
)



##################################################
# 8. SingleR annotation
##################################################

cat("Running SingleR annotation...\n")


ref <- HumanPrimaryCellAtlasData()


singleR_result <- SingleR(
    test=as.SingleCellExperiment(seu),
    ref=ref,
    labels=ref$label.main
)


seu$SingleR_annotation <- singleR_result$labels


write.csv(
    data.frame(
        cell=rownames(singleR_result),
        annotation=singleR_result$labels
    ),
    file.path(
        outdir,
        "SingleR_annotation.csv"
    ),
    row.names=FALSE
)



##################################################
# 9. Save plots
##################################################

p1 <- DimPlot(
    seu,
    reduction="umap",
    group.by="seurat_clusters",
    label=TRUE
)


ggsave(
    file.path(outdir,"UMAP_clusters.pdf"),
    p1,
    width=8,
    height=6
)



p2 <- DimPlot(
    seu,
    reduction="umap",
    group.by="SingleR_annotation",
    label=TRUE
)


ggsave(
    file.path(outdir,"UMAP_annotation.pdf"),
    p2,
    width=10,
    height=7
)



##################################################
# 10. Save Seurat object
##################################################

saveRDS(
    seu,
    file=file.path(
        outdir,
        "seurat_object.rds"
    )
)



cat("\nDONE\n")
