###################################################
#                                                 #
#   Convert Seurat object to annData that can be  #
#   accepted by GLUE                              #
#                                                 #
###################################################


# General libraries
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
  "SeuratDisk"
)
libraries(libs)


# Parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/scglue/"
# from.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/"
data.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scglue/"
setwd(data.dir)
getwd()
# args <- commandArgs(TRUE)
# file.path <- args[1]


# Conversion
message ("Began converting Seurat objects")
pbmclapply(read.table(paste0(code.dir, "../dataset_list.txt")) %>% `[[` (1), 
           mc.cores = detectCores(), function(file.path) {
  message ("----> Converting the file ", file.path, " ...")
  name <- strsplit(file.path, split = "/") %>% `[[` (1) %>% tail(., n = 1) %>% gsub(".qsave", "", .)
  name
  
  
  # Load data and generate scATAC-seq bed file
  obj <- qs::qread(file.path)
  dim(obj)
  
  
  # Convert the RNA and ATAC assays of Seurat object to Scanpy annData
  SaveH5Seurat(obj, filename = paste0(data.dir, name, ".h5Seurat"))
  Convert(paste0(data.dir, name, ".h5Seurat"), dest = "h5ad", assay = 'RNA',
          overwrite = TRUE)
  system(paste0("mv ", data.dir, name, ".h5ad", " ", data.dir, name,
                "-GEX.h5ad"))
  Convert(paste0(data.dir, name, ".h5Seurat"), dest = "h5ad", assay = 'ATAC',
          overwrite = TRUE)
  system(paste0("mv ", data.dir, name, ".h5ad", " ", data.dir, name,
                "-ATAC.h5ad"))
  
  
  # Convert the RNA and ATAC assays of Seurat object to csv files
  write.csv(obj[['RNA']]@counts, file = paste0(data.dir, name, "-GEX.csv"), 
            quote = FALSE)
  write.csv(obj[['ATAC']]@counts, file = paste0(data.dir, name, "-ATAC.csv"), 
            quote = FALSE)
})
message ("Finished converting Seurat objects")


save.image(paste0(data.dir, "1-GLUE-convert-rda-to-h5ad.RData"))
