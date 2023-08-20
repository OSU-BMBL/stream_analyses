# Load Packages
library(Seurat)
library(qs)
library(Seurat)
library(stringr)
library(Pando)
library(Seurat)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)

list <- readLines("List.txt")

# Get motif data
data('motifs')
data('motif2tf')

name <- list[1]
qs <- qread(file = paste0('data/', name,'.qsave'))

# Select variable features
seurat_object <- Seurat::FindVariableFeatures(qs, assay='RNA')

# extract gene annotations from EnsDb
#annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = FALSE))
#seqlevelsStyle(annotations) <- 'UCSC'
#Annotation(seurat_object@assays[["ATAC"]]) <- annotations

# Initiate GRN object and select candidate regions
seurat_object <- initiate_grn(seurat_object, rna_assay = 'RNA', peak_assay = 'ATAC')

# Scan candidate regions for TF binding motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

seurat_object <- find_motifs(
    seurat_object,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
)


# Infer gene regulatory network
seurat_object <- infer_grn(seurat_object)

# Print inferred coefficients
coef(seurat_object)

# Find gene and regulatory modules 
seurat_object <- find_modules(seurat_object)

# Print modules
modules <- NetworkModules(seurat_object) 
plot_gof(seurat_object, point_size=3)
plot_module_metrics(seurat_object)

## Visualizing the GRN
#seurat_object <- get_network_graph(seurat_object)
seurat_object <- get_network_graph(seurat_object,umap_method="corr")
qsave(seurat_object, file = paste0('Output/Pando/',name,'.qs'))

plot_network_graph(seurat_object)
plot_network_graph(seurat_object, layout='fr')