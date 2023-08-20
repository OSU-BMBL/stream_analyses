# Load Packages
library(Seurat)
library(qs)
library(Seurat)
library(stringr)
library(Pando)
library(Seurat)
library(tidyverse)
library(stats)

list <- readLines("List.txt")

# Get motif data

suppressMessages(library(ArchR))
suppressMessages(library(Signac))
suppressMessages(library(scMEGA))
suppressMessages(library(Nebulosa))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(EnsDb.Hsapiens.v75))
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(JASPAR2020))
suppressMessages(library(TFBSTools))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(MOJITOO))

# Get a list of motif position frequency matrices from the JASPAR database
reference <- readRDS("R/seurat.rds")
reference <- reference[, reference$celltype != "NA"]
reference <- reference %>%
    SCTransform(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
 Genome <- 'mm10'
name <- list[24]
qs <- qread(file = paste0('data/', name,'.qsave'))

qs@active.assay <- 'ATAC'
obj.atac <- qs
qs@active.assay <- 'RNA'
obj.rna <- qs

annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = FALSE))
seqlevelsStyle(annotations) <- 'UCSC'
# add the gene information to the object
Annotation(obj.atac) <- annotations

obj.rna <- obj.rna %>%
    SCTransform(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
obj.atac <- obj.atac %>%
    RunTFIDF() %>%
    FindTopFeatures() %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:30, verbose = FALSE)

# run sctransform

transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = obj.rna,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype,
  weight.reduction = obj.rna[['pca']],
  dims = 1:30
)

obj.rna <- AddMetaData(
  object = obj.rna,
  metadata = predictions
)
obj.atac <- AddMetaData(
  object = obj.atac,
  metadata = predictions
)

meta.data <- obj.rna@meta.data %>%
    as.data.frame()
pbmc <- CreateSeuratObject(
  counts = obj.rna@assays$RNA@counts,
  assay = "RNA",
    meta.data = meta.data
)
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = obj.atac@assays$ATAC@counts,
  sep = c(":", "-"),
    min.cells = 1,
    genome = Genome,
)
DefaultAssay(pbmc) <- "RNA"
pbmc <- pbmc %>%
    NormalizeData(verbose=F) %>%
    FindVariableFeatures(nfeatures=3000, verbose=F) %>%
    ScaleData(verbose=F) %>%
    RunPCA(npcs=50, reduction.name="RNA_PCA", verbose=F)
DefaultAssay(pbmc) <- "ATAC"
pbmc <- pbmc %>%
    RunTFIDF(verbose=F) %>%
    FindTopFeatures(min.cutoff = 'q0', verbose=F) %>%
    RunSVD(verbose=F)

pbmc <- mojitoo(
     object = pbmc,
     reduction.list = list("RNA_PCA", "lsi"),
     dims.list = list(1:50, 2:50), ## exclude 1st dimension of LSI
     reduction.name = 'MOJITOO',
     assay = "RNA"
)
DefaultAssay(pbmc) <- "RNA"
embedd <- Embeddings(pbmc[["MOJITOO"]])
pbmc <- RunUMAP(pbmc, 
                reduction="MOJITOO", 
                reduction.name="MOJITOO_UMAP", 
                dims=1:ncol(embedd), verbose=F)

D1 <- DimPlot(pbmc, group.by = "predicted.id", 
        shuffle = TRUE, label = TRUE, reduction = "MOJITOO_UMAP") + NoLegend()
# Trajectory analysis
pbmc <- AddTrajectory(object = pbmc, 
                      trajectory = c("naive CD4 T cells", 
                                     "memory B cells"),
                      group.by = "predicted.id", 
                          reduction = "MOJITOO_UMAP",
                          dims = 1:2, 
                          use.all = FALSE)
pbmc.t.cells <- pbmc[, !is.na(pbmc$Trajectory)]

p1 <- DimPlot(object = pbmc.t.cells, 
              group.by = "predicted.id", 
              reduction = "MOJITOO_UMAP",
             label = TRUE) + NoLegend()

p2 <- TrajectoryPlot(object = pbmc.t.cells, 
                    reduction = "MOJITOO_UMAP",
                    continuousSet = "blueYellow",
                    size = 1,
                   addArrow = FALSE)
# Select TFs
pbmc.t.cells <- AddMotifs(
  object = pbmc.t.cells,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm,
    assay = "ATAC")
pbmc.t.cells <- RunChromVAR(
  object = pbmc.t.cells,
  genome = BSgenome.Mmusculus.UCSC.mm10,
    assay = "ATAC")

sel.tfs <- SelectTFs(object = pbmc.t.cells)
df.cor <- sel.tfs$tfs
ht <- sel.tfs$heatmap
draw(ht)

# Select genes
sel.genes <- SelectGenes(object = pbmc.t.cells,
                  labelTop1 = 0,
                  labelTop2 = 0)
df.p2g <- sel.genes$p2g
ht <- sel.genes$heatmap
draw(ht)

# Gene regulatory network inference and visualization
tf.gene.cor <- GetTFGeneCorrelation(object = pbmc.t.cells, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")
ht <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point)
ht
time_point <- df.cor$time_point
save(list=c("tf.gene.cor", "ht", "time_point"), file=paste0('STREAM Revision/Output/scMEGA/',name,'.RData'))
