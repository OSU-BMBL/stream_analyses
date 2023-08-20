######################################################
#                                                    #
#   Preprocess the scATAC + scGEX data set from      #
#   SNARE-seq2 on human GM12878 and A549 samples     #
#                                                    #
######################################################


# Libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(easypackages)
libs <- c(
  "qs", 
  "hdf5r",
  "Seurat",
  "Signac", 
  "pbmcapply", 
  "pbapply", 
  "parallel",
  "EnsDb.Mmusculus.v79", 
  "EnsDb.Hsapiens.v86",
  "dplyr", 
  "ggplot2", 
  "ggpubr", 
  "dplyr"
)
libraries(libs)


# Parameters
code.dir <- "/fs/ess/PCON0022/liyang/Joint-ATAC-RNA/SNARE-seq2/"
data.dir <- "/fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/SNARE-seq2/rda/"
dir.create("/fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/SNARE-seq2/rda/")
org <- "hg38"
setwd(data.dir)
getwd()


# load the RNA data
rna_counts <- Read10X(data.dir = "../human.snare2-rna.matrix")
dim(rna_counts)
rna_counts[1:4, 1:4]


# load the ATAC data
atac_counts <- Read10X(data.dir = "../human.snare2-atac.matrix")
dim(atac_counts)
atac_counts[1:4, 1:4]


# Unify cells between GEX and ATAC
common.cells <- intersect(colnames(rna_counts), colnames(atac_counts))
rna_counts <- rna_counts[, common.cells]
atac_counts <- atac_counts[, common.cells]
dim(rna_counts)
dim(atac_counts)
identical(colnames(rna_counts), colnames(atac_counts))


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


# create a Seurat object containing the RNA data
pbmc <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA"
)
dim(pbmc)


# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  annotation = annotation
)
dim(pbmc[["ATAC"]])


# split the dataset into a list of two seurat objects (GM12878 and A549)
ifnb.list <- SplitObject(pbmc, split.by = "orig.ident")
names(ifnb.list)
A549 <- ifnb.list[[1]]
GM <- ifnb.list[[2]]


# Quality control
# DefaultAssay(pbmc) <- "ATAC"
# pbmc <- NucleosomeSignal(pbmc)
# pbmc <- TSSEnrichment(pbmc)
DefaultAssay(A549) <- "RNA"
A549[["percent.mt"]] <- PercentageFeatureSet(A549, pattern = "^MT-")
qs::qsave(A549, paste0(data.dir, "A549_Seu_before_QC.qsave"))

DefaultAssay(GM) <- "RNA"
GM[["percent.mt"]] <- PercentageFeatureSet(GM, pattern = "^MT-")
qs::qsave(GM, paste0(data.dir, "GM_Seu_before_QC.qsave"))


# # Read metadata
# meta.df <- read.csv("../GSE214979_cell_metadata.csv", row.names = 1)
# head(meta.df)
# dim(meta.df)
# identical(meta.df, pbmc@meta.data)


png(paste0(data.dir, "A549_RNA_ATAC_MT_before_QC.png"), 
    res = 150, 
    width = 2000, height = 1000)
VlnPlot(
  object = A549,
  features = c("nCount_RNA", "nCount_ATAC", 
               "percent.mt"),
  ncol = 4,
  pt.size = 0
)
dev.off()

png(paste0(data.dir, "GM_RNA_ATAC_MT_before_QC.png"), 
    res = 150, 
    width = 2000, height = 1000)
VlnPlot(
  object = GM,
  features = c("nCount_RNA", "nCount_ATAC", 
               "percent.mt"),
  ncol = 4,
  pt.size = 0
)
dev.off()


# filter out low quality cells
A549 <- subset(
  x = A549,
  subset = nCount_RNA < 15000 &
    nCount_RNA > 200 & 
    nCount_ATAC < 3500 & 
    nCount_ATAC > 200 & 
    # nucleosome_signal < 1 &
    # TSS.enrichment > 2 & 
    percent.mt < 5
)
dim(A549)
qs::qsave(A549, paste0(data.dir, "SNAREseq2_hg38_A549.qsave"))
system("cp SNAREseq2_hg38_A549.qsave ..")


GM <- subset(
  x = GM,
  subset = nCount_RNA < 30000 &
    nCount_RNA > 200 & 
    nCount_ATAC < 10000 & 
    nCount_ATAC > 200 & 
    # nucleosome_signal < 1 &
    # TSS.enrichment > 2 & 
    percent.mt < 2
)
dim(GM)
qs::qsave(GM, paste0(data.dir, "SNAREseq2_hg38_GM12878.qsave"))
system("cp SNAREseq2_hg38_GM12878.qsave ..")


save.image(paste0(data.dir, "	1-Preprocess-SNARE-seq2-human-GM12878-A549.RData"))
