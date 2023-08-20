################################################################
#                                                              #
#                      TF-region calculation                   #
#                                                              #
################################################################


# General libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(stream2)
message ("Loaded libraries: ", paste(libs, collapse = ", "))


# Set parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/"
parent.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/"
eval.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/"
work.dir <- parent.dir
setwd(work.dir)
getwd()


# Get data set list
data.ls <- read.table(paste0(code.dir, "dataset_names.txt")) %>% `[[` (1)
data.ls
length(data.ls)
message ("There are ", length(data.ls), " datasets")


################################################################
#                                                              #
#                   Get TF-region for six methods              #
#                                                              #
################################################################


#
# STREAM
#
message ("Began converting eRegulons from nested lists to data.table")
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/stream/"
setwd(work.dir)
getwd()
stream.tfRegion <- pbmclapply(seq_along(data.ls), mc.cores = detectCores(), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.qsave"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res <- qs::qread(paste0(work.dir, data.ls[i], "/eRegulons.qsave")) %>% lapply(., function(x) {
    data.frame(
      TF = rep(x$TF, length(x$peaks)),
      Region = x$peaks
    )
  }) %>% Reduce("rbind", .)
  return (res[!duplicated(res), , drop = FALSE])
})
head(stream.tfRegion[[1]])
names(stream.tfRegion) <- data.ls
qs::qsave(stream.tfRegion, paste0(parent.dir, "stream_TF_region_links.qsave"))
setwd("..")


#
# SCENIC+
#
# Conversion
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scenicplus/"
setwd(work.dir)
scplus.tfRegion <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  dt <- read.csv(paste0(work.dir, data.ls[i], "/eRegulons.csv"), sep = ",")[, c("TF", "Region")]
  return(dt[!duplicated(dt),, drop = FALSE])
})
head(scplus.tfRegion)
names(scplus.tfRegion) <- data.ls
scplus.tfRegion <- scplus.tfRegion[!sapply(scplus.tfRegion, is.null)]
length(scplus.tfRegion)
qs::qsave(scplus.tfRegion, paste0(parent.dir, "scenicplus_TF_region_links.qsave"))
setwd("..")


#
# Pando
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/pando/"
setwd(work.dir)
library(Pando)
pando.tfRegion <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".qs"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn <- qs::qread(paste0(work.dir, data.ls[i], ".qs")) %>% NetworkModules
  dt <- grn@meta[, c("tf", "regions")]
  colnames(dt) <- c("TF", "Region")
  return(dt[!duplicated(dt),, drop = FALSE])
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(pando.tfRegion[[1]])
names(pando.tfRegion) <- data.ls
pando.tfRegion <- pando.tfRegion[!sapply(pando.tfRegion, is.null)]
length(pando.tfRegion)
qs::qsave(pando.tfRegion, paste0(parent.dir, "pando_TF_region_links.qsave"))
setwd("..")


#
# GLUE
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scglue/"
setwd(work.dir)
scglue.tfRegion <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn <- read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  colnames(grn) <- c("TF", "Gene", "importance")
  gene2peak <- read.csv(paste0(work.dir, data.ls[i], "/gene2peak.csv"), sep = ",")
  colnames(gene2peak) <- c("Gene", "Region", "Weight")
  dt <- merge(grn, gene2peak, by = "Gene")
  dt <- dt[, c("TF", "Region", "Gene")]
  return(dt[!duplicated(dt),, drop = FALSE])
})
head(scglue.tfRegion[[1]])
names(scglue.tfRegion) <- data.ls
scglue.tfRegion <- scglue.tfRegion[!sapply(scglue.tfRegion, is.null)]
length(scglue.tfRegion)
qs::qsave(scglue.tfRegion, paste0(parent.dir, "scglue_TF_region_links.qsave"))
setwd("..")


# Get list of predicted TF-to-region linkages for each method
method.list <- list(
  stream = stream.tfRegion, 
  scenicplus = scplus.tfRegion, 
  scglue = scglue.tfRegion,
  pando = pando.tfRegion
)


################################################################
#                                                              #
#        Benchmark TF-region linkages for each data set        #
#                                                              #
################################################################


# Read correspondence file
top.peaks <- 10000
chipseq.hic <- read.csv(file = paste0(parent.dir, "ChIPseq-and-HiC-data.csv"))
dim(chipseq.hic)
head(chipseq.hic)


# Read TF lists from each TF ChIP-seq file
chipseq.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/"
message ("Begin reading TF ChIP-seq file names")
pbmclapply(1:nrow(chipseq.hic), mc.cores = detectCores(), function(i) {
  message ("----> Reading TF ChIP-seq for the ", i, "-th dataset: ", chipseq.hic$Dataset[i])
  if (chipseq.hic$ChIP.seq[i] == "N/A") {
    return(NULL)
  }
  dir.ll <- strsplit(chipseq.hic$ChIP.seq[i], split = ", ") %>% `[[` (1)
  lapply(dir.ll, function(x) {
    system( paste0("gunzip ", eval.dir, x, "/*.bed.gz") )
    message ("--------> unziped .bed.gz files in dir: ", x)
  })
  res <- lapply(dir.ll, function(x) {
    message("Reading bed files in dir: ", x)
    bed.files <- list.files(paste0(eval.dir, x), pattern = ".bed")
    dir.tfRegion <- pblapply(bed.files, function(y) {
      bed.dt <- read.table(paste0(eval.dir, x, "/", y)) %>% arrange(desc(V7)) %>% dplyr::select(V1, V2, V3) %>% 
        head(n = top.peaks)
      if (grepl("hg", x)) {
        TF <- toupper(strsplit(y, split = "_")[[1]][1])
      }
      data.frame(
        TF = rep(TF, nrow(bed.dt)),
        region = paste0(bed.dt$V1, "-", bed.dt$V2, "-", bed.dt$V3)
      )
    }) %>% Reduce("rbind", .)
  }) %>% Reduce("rbind", .)
  qs::qsave(res, paste0(eval.dir, chipseq.hic$Dataset[i], "_ChIPseq.qsave"))
})


################################################################
#                                                              #
#                     TF-to-region calculation                 #
#                                                              #
################################################################


# Steps:
# 1. Select the datasets composed of at least ten TFs
# 2. For each TF in each dataset, calculate precision, recall, and f scores
# 3. For each dataset, draw box plots for precision, recall, and f scores
# 4. Compare benchmarking methods


# Select data containing more than ten TFs
#tf.covData <- qs::qread(paste0(parent.dir, "data-list-containing-more-than-10-tfs.qsave"))
tf.covData <- chipseq.hic$Dataset
tf.covData


# Define the function to calculate region overlaps
cal_overlap_regions <- function(r1, r2) {
  
  if (is.null(r1) | is.null(r2)) {
    return(c(precision = 0, recall = 0, f.score = 0))
  }
  invisible(library(GenomicRanges))
  invisible(library(Signac))
  hits <- findOverlaps(StringToGRanges(r1, sep = c(":", "-")), 
                       StringToGRanges(r2, sep = c(":", "-")) )
  prec <- length(unique(queryHits(hits))) / length(r1)
  recall <- length(unique(subjectHits(hits))) / length(r2)
  c(precision = prec, recall = recall, f.score = 2 * prec * recall / (prec + recall))
}


# Define function to calculate precision, recall, and f-score for TF-to-region prediction
# x : the prediction vector
# y : the labels 
cal_prec_recall_fscore <- function(x, y) {
  
  tfs <- intersect(unique(x$TF), unique(y$TF) )
  message ("--------> ", length(tfs), " common TFs to calculate TF-to-region linkages")
  if (length(tfs) < 1) {
    return(data.frame(precision = 0, recall = 0, f.score = 0))
  }
  tf.scores <- lapply(tfs, function(z) {
    x.regions <- unique(x[x$TF == z,, drop = FALSE]$Region) %>% strsplit(., split = ";") %>% 
      sapply(., "[[", 1)
    y.regions <- unique(y[y$TF == z,, drop = FALSE]$region)
    cal_overlap_regions( r1 = x.regions, r2 = y.regions )
  }) %>% Reduce("rbind", .)
  if (is.vector(tf.scores)) {
    tf.scores <- t(as.data.frame(tf.scores) )
  }
  rownames(tf.scores) <- tfs
  as.data.frame(tf.scores)
}


# Calculate three scores of TF-to-region scores for all methods
library(data.table)
res.ll <- pbmclapply(tf.covData, mc.cores = detectCores(), function(chip.data) {
  message ("----> ChIP-seq data: ", chip.data)
  lapply(seq_along(method.list), function(i) {
    message (i)
    method <- method.list[[i]]
    if (!chip.data %in% names(method) | !file.exists(paste0(chipseq.dir, chip.data, "_ChIPseq.qsave"))) {
      return(c(prec = 0, recall = 0, f.score = 0, method = names(method.list)[i]))
    }
    method.scores <- cal_prec_recall_fscore(x = method[[chip.data]], 
                                            y = qs::qread(paste0(chipseq.dir, chip.data, "_ChIPseq.qsave")) )
    method.scores[is.na(method.scores)] <- 0.0
    cbind(method.scores, method = rep(names(method.list)[i], nrow(method.scores)))
  }) %>% Reduce("rbind", .)
})
sapply(res.ll, dim)
names(res.ll) <- tf.covData
qs::qsave(res.ll, paste0(parent.dir, "TF-to-region-benchmark-all.qsave"))


# Select representative datasets
mean.prec <- sapply(seq_along(res.ll), function(i) {
  message (i)
  x <- res.ll[[i]]
  split(x, f = x$method) %>% sapply(., function(xx) {
    mean(as.numeric(xx$precision) )
  })
})
colnames(mean.prec) <- names(res.ll)
selected.data <- names(res.ll)
selected.scores <- mean.prec[, selected.data]
intersect(qs::qread(paste0(parent.dir, "TF-coverage-benchmark.qsave")) %>% names, 
          selected.data)
qs::qsave(selected.scores, paste0(parent.dir, "TF-to-region-benchmark.qsave"))


save.image("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/2-benchmark-TF-region.RData")
# load("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/2-benchmark-TF-region.RData")