######################################################
#                                                    #
#   Preprocess the scATAC + scGEX data set from      #
#   10X Genomics Multiome on human bone marrow       #
#                                                    #
######################################################


# # Environment variables
# PATH <- "/users/PAS1475/liyang/libs/hdf5_1.10.6/bin:/apps/xalt/xalt/bin:/opt/mvapich2/intel/19.0/2.3.3/bin:/apps/gnu/8.4.0/bin:/opt/intel/19.0.5/itac_latest/bin:/opt/intel/19.0.5/advisor/bin64:/opt/intel/19.0.5/vtune_amplifier/bin64:/opt/intel/19.0.5/inspector_2019/bin64:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/bin/intel64:/usr/lib64/qt-3.3/bin:/opt/osc/bin:/usr/local/bin:/usr/bin:/app/cmake/3.25.2/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/fs/ess/PCON0022/liyang/tools/MEME/bin:/fs/ess/PCON0022/liyang/tools/MEME/libexec/meme-5.5.0:/fs/ess/PCON0022/liyang/tools/HOMER/bin/:/apps/cmake/3.25.2/bin/"
# LD_LIBRARY_PATH <- "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib:/opt/mvapich2/intel/19.0/2.3.3/lib:/apps/gnu/8.4.0/lib64:/apps/gnu/8.4.0/lib:/opt/intel/19.0.5/debugger_2019/libipt/intel64/lib:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/lib/intel64_lin:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/daal/lib/intel64_lin:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/ipp/lib/intel64_lin:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/mkl/lib/intel64_lin:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/tbb/lib/intel64_lin/gcc4.4:/opt/ddn/cci/lib:/opt/ddn/ime/lib:/opt/ddn/isa-l/lib"
# CPATH <- "/users/PAS1475/liyang/libs/hdf5_1.10.6/include:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/mkl/include:/opt/intel/19.0.5/compilers_and_libraries_2019/linux/tbb/include"
# Sys.setenv(PATH = PATH)
# Sys.setenv(LD_LIBRARY_PATH = LD_LIBRARY_PATH)
# Sys.setenv(CPATH = CPATH)
# Sys.getenv("PATH")
# Sys.getenv("LD_LIBRARY_PATH")
# Sys.getenv("CPATH")


# Libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(easypackages)
libs <- c(
  "qs", 
  "hdf5r",
  "Seurat",
  "SeuratDisk",
  "Signac", 
  "pbmcapply", 
  "pbapply", 
  "parallel",
  "EnsDb.Mmusculus.v79", 
  "EnsDb.Hsapiens.v86",
  "dplyr", 
  "ggplot2", 
  "ggpubr"
)
libraries(libs)


# Parameters
code.dir <- "/fs/ess/PCON0022/liyang/Joint-ATAC-RNA/10X-human-bone-marrow/"
data.dir <- "/fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/10X-Genomics-Multiome/10X-human-bone-marrow/rda/"
dir.create("/fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/10X-Genomics-Multiome/10X-human-bone-marrow/rda/")
org <- "hg38"
setwd(data.dir)
getwd()


# Convert h5ad to h5Seurat for GEX assay
Convert("../bm_multiome_rna.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc_rna <- LoadH5Seurat("../bm_multiome_rna.h5seurat")
str(pbmc_rna)
rownames(pbmc_rna)
colnames(pbmc_rna)


# Convert gene names from integers to symbols
rownames(pbmc_rna[['RNA']]@counts) <- read.table(paste0(data.dir, "../bm_multiome_rna.txt"))[[1]]
rownames(pbmc_rna[['RNA']]@data) <- rownames(pbmc_rna[['RNA']]@counts)
rownames(pbmc_rna[['RNA']]@meta.features) <- rownames(pbmc_rna[['RNA']]@counts)


# Convert h5ad to h5Seurat for ATAC assay
Convert("../bm_multiome_atac.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc_atac <- LoadH5Seurat("../bm_multiome_atac.h5seurat")
str(pbmc_atac)
rownames(pbmc_atac)
colnames(pbmc_atac)


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


# Construct Seurat object composed of both ATAC and GEX
pbmc <- pbmc_rna
pbmc[['ATAC']] <- CreateChromatinAssay(
    counts = pbmc_atac@assays$RNA@counts,
    sep = c(":", "-"),
    annotation = annotation
  )
str(pbmc)
head(pbmc@meta.data)
dim(pbmc@meta.data)
colnames(pbmc@meta.data)
pbmc@meta.data[1:4, 1:4]
pbmc[['RNA']][1:4, 1:4]
pbmc[['ATAC']][1:4, 1:4]


qs::qsave(pbmc, paste0(data.dir, "10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave"))
system("cp 10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave ../..")
save.image(paste0(data.dir, "2-Preprocess-10X-human-bone-marrow.Rdata"))
# load("/fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/10X-Genomics-Multiome/10X-human-bone-marrow/rda/2-Preprocess-10X-human-bone-marrow.Rdata")