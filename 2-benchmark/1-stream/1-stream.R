###################################################
#                                                 #
#       Run STREAM on benchmarking dataset        #
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
  "igraph", 
  "sparseMatrix"
)
libraries(libs)


# Specialized libraries
# library(stream2)
source.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream_v2/R/"
code.list <- list.files(source.dir, pattern = ".R")
lapply(code.list, function(x) {
  message (x)
  source(paste0(source.dir, x))
})


# Parameters
args <- commandArgs(TRUE)
file.path <- args[1]
file.name <- strsplit(file.path, split = "\\/") %>% `[[` (1) %>% tail(n = 1)
data.dir <- paste0("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/stream/", 
                   gsub(".qsave", "", file.name))
dir.create(data.dir)
data.dir
org <- strsplit(file.name, split = "_") %>% `[[` (1) %>% `[` (2)
org


# Run STREAM
setwd(data.dir)
getwd()
obj <- qs::qread(file.path)
dim(obj)
en.regs <- run_stream(obj = obj, org = org)
save.image(paste0(data.dir, "/", gsub(".qsave", "", file.name), ".RData"))
qs::qsave(en.regs, paste0(data.dir, "/eRegulons.qsave"))