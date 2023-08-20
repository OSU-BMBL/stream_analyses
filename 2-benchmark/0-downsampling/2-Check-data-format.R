#####################################################
#                                                   #
#           Check the formats of datasets           #
#                                                   #
#####################################################


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
  "stringr"
)
libraries(libs)


# Parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Downsampling/"
data.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/"
setwd(data.dir)
getwd()


# 10X Genomics Multiome ATAC+GEX
setwd("10X/")
dir(pattern = ".qsave")


# DLPFC
x <- qs::qread("10X_hg38_DLPFC.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Human healthy brain
x <- qs::qread("10X_hg38_Flash_Frozen_Human_Healthy_Brain_Tissue_3k.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Human bone marrow
x <- qs::qread("10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Mouse brain
x <- qs::qread("10X_mm10_Mouse_brain_nuclei_isolated.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Mouse kidney
x <- qs::qread("10X_mm10_Mouse_kidney_nuclei_isolated.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# ASTAR-seq
setwd("../ASTARseq")
getwd()
dir(pattern = ".qsave")


# K562
x <- qs::qread("ASTAR_seq_hg19_K562.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# K562 + BJ
x <- qs::qread("ASTARSeq_hg19_K562_BJ.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# NEAT-seq
setwd("../NEATseq/")
dir(pattern = ".qsave")


# GM12878
x <- qs::qread("NEATseq_hg38_GM12878.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# K562
x <- qs::qread("NEATseq_hg38_K562.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# mESCs
x <- qs::qread("NEATseq_mm10_mESCs.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# scCAT-seq
setwd("../scCATseq/")
dir(pattern = ".qsave")


# HCT116 + K562 + HelaS3 + PDX
x <- qs::qread("scCAT_seq_hg19_K562_HCT116_HelaS3_PDX.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Implanted embryos
x <- qs::qread("scCAT_seq_hg19_PreImplanted_embryos.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# sciCAR
setwd("../sciCAR/")
dir(pattern = ".qsave")


# A549
x <- qs::qread("sciCAR_hg19_A549.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# HEK293T
x <- qs::qread("sciCAR_hg19_HEK293T.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Mouse healthy kidney
x <- qs::qread("sciCAR_mm10_Mouse_healthy_kidney.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# scNMT-seq
setwd("../scNMTseq")
dir(pattern = ".qsave")


# mESCs
x <- qs::qread("scNMT_seq_mm10_mESCs.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Mouse embryos
x <- qs::qread("scNMT_seq_mm10_Mouse_embryos.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# SHARE-seq
setwd("../SHAREseq")
dir(pattern = ".qsave")


# GM12878 mixed with NIH/3T3
x <- qs::qread("SHAREseq_hg19_GM12878_mixed_with_NIH3T3.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# GM12878 Rep1
x <- qs::qread("SHAREseq_hg19_GM12878_rep1.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# GM12878 Rep2
x <- qs::qread("SHAREseq_hg19_GM12878_rep2.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# GM12878 Rep3
x <- qs::qread("SHAREseq_hg19_GM12878_rep3.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# K562 mixed with mouse RAW cells
x <- qs::qread("SHAREseq_hg19_K562_mixed_with_mouse_RAW_cells.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Adult mouse lung
x <- qs::qread("SHAREseq_mm10_Adult_mouse_lung.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# Mouse brain
x <- qs::qread("SHAREseq_mm10_Mouse_brain.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# SNARE-seq2
setwd("../SNAREseq2")
dir(pattern = ".qsave")


# A549
x <- qs::qread("SNAREseq2_hg38_A549.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")


# GM12878
x <- qs::qread("SNAREseq2_hg38_GM12878.qsave")
x[['RNA']][1:3, 1:3]
x[['ATAC']][1:3, 1:3]
dir(pattern = ".qsave")
