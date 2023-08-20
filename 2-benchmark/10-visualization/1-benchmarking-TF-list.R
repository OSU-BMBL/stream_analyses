################################################################
#                                                              #
#                    Convert results to RData                  #
#                                                              #
################################################################


# General libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(stream2)
message ("Loaded libraries: ", paste(libs, collapse = ", "))


# Set parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/"
parent.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/"
work.dir <- parent.dir
setwd(work.dir)
getwd()


# Get scenicplus results
data.ls <- read.table(paste0(code.dir, "dataset_names.txt")) %>% `[[` (1)
data.ls
length(data.ls)
message ("There are ", length(data.ls), " datasets")


################################################################
#                                                              #
#                   Get TF list for six methods                #
#                                                              #
################################################################


# 
# STREAM
#
message ("Began converting eRegulons from nested lists to data.table")
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/stream/"
setwd(work.dir)
getwd()
stream.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.qsave"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  return (qs::qread(paste0(work.dir, data.ls[i], "/eRegulons.qsave")) %>% lapply(., "[[", "TF") %>% 
    unlist %>% unique)
  message ("Number of eRegulons: ", length(regs))
  
})
head(stream.tfs)
names(stream.tfs) <- data.ls
qs::qsave(stream.tfs, paste0(parent.dir, "stream_TFs.qsave"))
setwd("..")


#
# SCENIC+
#
# Conversion
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scenicplus/"
setwd(work.dir)
scplus.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  dt <- read.csv(paste0(work.dir, data.ls[i], "/eRegulons.csv"), sep = ",")
  return(unique(split(dt, f = dt$TF) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
})
head(scplus.tfs)
names(scplus.tfs) <- data.ls
scplus.tfs <- scplus.tfs[!sapply(scplus.tfs, is.null)]
length(scplus.tfs)
qs::qsave(scplus.tfs, paste0(parent.dir, "scenicplus_TFs.qsave"))
setwd("..")


#
# SCENIC
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scenic/"
setwd(work.dir)
scenic.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  dt <- read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  return(unique(split(dt, f = dt$TF) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(scenic.tfs)
names(scenic.tfs) <- data.ls
scenic.tfs <- scenic.tfs[!sapply(scenic.tfs, is.null)]
length(scenic.tfs)
qs::qsave(scenic.tfs, paste0(parent.dir, "scenic_TFs.qsave"))
setwd("..")


# 
# DIRECT-NET
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/directnet/"
setwd(work.dir)
dinet.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/", data.ls[i], 
                          ".qsave_directnet_out.qs"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  df <- qs::qread(paste0(work.dir, data.ls[i], "/", data.ls[i], 
                        ".qsave_directnet_out.qs"))
  dt <- df$link
  return(unique(split(dt, f = dt$TF) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(dinet.tfs)
names(dinet.tfs) <- data.ls
dinet.tfs <- dinet.tfs[!sapply(dinet.tfs, is.null)]
length(dinet.tfs)
qs::qsave(dinet.tfs, paste0(parent.dir, "directnet_TFs.qsave"))
setwd("..")


#
# Pando
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/pando/"
setwd(work.dir)
library(Pando)
pando.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".qs"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn <- qs::qread(paste0(work.dir, data.ls[i], ".qs")) %>% NetworkModules
  dt <- grn@meta
  return(unique(split(dt, f = dt$tf) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(pando.tfs)
names(pando.tfs) <- data.ls
pando.tfs <- pando.tfs[!sapply(pando.tfs, is.null)]
length(pando.tfs)
qs::qsave(pando.tfs, paste0(parent.dir, "pando_TFs.qsave"))
setwd("..")


#
# GLUE
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scglue/"
setwd(work.dir)
scglue.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn <- read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  gene2peak <- read.csv(paste0(work.dir, data.ls[i], "/gene2peak.csv"), sep = ",")
  grn <- grn[grn$target %in% unique(gene2peak$Source),, drop = FALSE]
  dt <- grn
  return(unique(split(dt, f = dt$TF) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(scglue.tfs)
names(scglue.tfs) <- data.ls
scglue.tfs <- scglue.tfs[!sapply(scglue.tfs, is.null)]
length(scglue.tfs)
qs::qsave(scglue.tfs, paste0(parent.dir, "scglue_TFs.qsave"))
setwd("..")


#
# scMEGA
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scMEGA/"
setwd(work.dir)
scmega.tfs <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  dt <- read.csv(paste0(work.dir, data.ls[i], ".csv"), sep = ",")
  return(unique(split(dt, f = dt$tf) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(scmega.tfs)
names(scmega.tfs) <- data.ls
scmega.tfs <- scmega.tfs[!sapply(scmega.tfs, is.null)]
length(scmega.tfs)
qs::qsave(scmega.tfs, paste0(parent.dir, "scmega_TFs.qsave"))
setwd("..")


################################################################
#                                                              #
#                   Benchmark TF for each data set             #
#                                                              #
################################################################


# Read correspondence file
chipseq.hic <- read.csv(file = paste0(parent.dir, "ChIPseq-and-HiC-data.csv"))
dim(chipseq.hic)
head(chipseq.hic)


# Read TF lists from each TF ChIP-seq file
chipseq.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/"
message ("Begin reading TF ChIP-seq file names")
chipseq.tfs <- pblapply(1:nrow(chipseq.hic), function(i) {
  message ("----> Reading TF ChIP-seq for the ", i, "-th dataset: ", chipseq.hic$Dataset[i])
  if (chipseq.hic$ChIP.seq[i] == "N/A") {
    return(NULL)
  }
  if (!grepl(",", chipseq.hic$ChIP.seq[i]) ) {
    tf.lst <- list.files(path = paste0(chipseq.dir, chipseq.hic$ChIP.seq[i]) ) %>% strsplit(., split = "_") %>% 
      sapply(., "[", 1) %>% unique
  } else {
    tf.lst <- sapply(strsplit(chipseq.hic$ChIP.seq[i], split = ", ") %>% `[[` (1), function(x) {
      list.files(path = paste0(chipseq.dir, x) ) %>% strsplit(., split = "_") %>% 
        sapply(., "[", 1) %>% unique
    }) %>% Reduce("union", .)
  }
  if (grepl("hg", chipseq.hic$Dataset[i])) {
    message ("----> Converting TF symbols to upper case")
    tf.lst <- toupper(tf.lst)
  }
  tf.lst
})
names(chipseq.tfs) <- chipseq.hic$Dataset
head(chipseq.tfs)
length(stream.tfs)
chipseq.tfs <- chipseq.tfs[!sapply(chipseq.tfs, is.null)]
length(chipseq.tfs)
chipseq.tfs <- chipseq.tfs[sapply(chipseq.tfs, length) > 0]
length(chipseq.tfs)
sapply(chipseq.tfs, length)
qs::qsave(chipseq.tfs, paste0(parent.dir, "ENCODE-TF-list.qsave"))
message ("Finish reading TF ChIP-seq file names")


################################################################
#                                                              #
#                      TF coverage calculation                 #
#                                                              #
################################################################


# Steps:
# 1. Select the datasets composed of at least ten TFs
# 2. Rank TFs in decreasing order of eRegulon/regulon ranks or degrees
# 3. Select the top-N TFs to ensure all methods could deliver the same numbers of TFs
# 4. For each dataset, calculate the number of covered TFs in TF ChIP-seq datasets; draw a curve
# 5. Calculate AUROC for each method


# Select data containing more than ten TFs
sapply(chipseq.tfs, length) %>% range
tf.covData <- names(chipseq.tfs)
tf.covData
qs::qsave(tf.covData, paste0(parent.dir, "data-list.qsave"))


# Select the top-N TFs to ensure all methods could deliver the same numbers of TFs
sapply(chipseq.tfs, length)
sapply(chipseq.tfs[tf.covData], length)
stream.tfs[intersect(tf.covData, names(stream.tfs))]
scplus.tfs[intersect(tf.covData, names(scplus.tfs))]
scenic.tfs[intersect(tf.covData, names(scenic.tfs))]
scglue.tfs[intersect(tf.covData, names(scglue.tfs))]
pando.tfs[intersect(tf.covData, names(pando.tfs))]
scmega.tfs[intersect(tf.covData, names(scmega.tfs))]
dinet.tfs[intersect(tf.covData, names(dinet.tfs))]
method.list <- list(
  stream = stream.tfs, 
  scenicplus = scplus.tfs, 
  scenic = scenic.tfs,
  scglue = scglue.tfs,
  pando = pando.tfs,
  scmega = scmega.tfs,
  directnet = dinet.tfs
)


# Define function to calculate precision, recall, and f-score for TF coverage
# x : the prediction vector
# y : the labels 
cal_prec_recall_fscore <- function(x, y) {
  tp <- intersect(x, y)
  prec <- length(tp) / length(x)
  recall <- length(tp) / length(y)
  c(
    prec = prec, 
    recall = recall, 
    f.score = 2 * prec * recall / (prec + recall)
  )
}


# Calculate three scores of TF coverage for all methods
library(data.table)
res.ll <- pblapply(tf.covData, function(chip.data) {
  message ("----> ChIP-seq data: ", chip.data)
  res <- lapply(method.list, function(method) {
    if (!chip.data %in% names(method)) {
      return(c(prec = 0, recall = 0, f.score = 0))
    }
    cal_prec_recall_fscore(method[[chip.data]], chipseq.tfs[[chip.data]])
  }) %>% Reduce("rbind", .)
  rownames(res) <- names(method.list)
  res[is.na(res)] <- 0
  res
})
names(res.ll) <- tf.covData


# Select representative datasets
benchTF.data <- tf.covData


save.image("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/1-benchmarking-TF-list.RData")
# load("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/1-benchmarking-TF-list.RData")