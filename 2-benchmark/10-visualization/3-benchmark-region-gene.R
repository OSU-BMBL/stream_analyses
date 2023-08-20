################################################################
#                                                              #
#                      Region-gene calculation                 #
#                                                              #
################################################################


# General libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(stream2)
message ("Loaded libraries: ", paste(libs, collapse = ", "))


# Set parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/"
parent.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/"
eval.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/Chrom-interaction-validation/"
gtf.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Feasibility/gtf-files/"
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
#                   Get region-gene for six methods            #
#                                                              #
################################################################


# Collapse transcripts of a gene to the longest one
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}


# Define the function to get gene locations
get_gene_loc <- function(symbols, org = "hg38") {
  
  # Load symbol-to-loci correspondence for genes
  if (org == "mm10") {
    library(EnsDb.Mmusculus.v79)
    annot <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79) %>% CollapseToLongestTranscript
  } else if (org == "hg19") {
    library(EnsDb.Hsapiens.v75)
    annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v75) %>% CollapseToLongestTranscript
  } else {
    library(EnsDb.Hsapiens.v86)
    annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86) %>% CollapseToLongestTranscript
  }
  message ("Loaded loci for genes of ", org)
  
  
  # Get locations for the input gene symbols
  invisible(library(Signac) )
  named.loci <- setNames(
    paste0("chr", GRangesToString(annot)), 
    annot$gene_name
  )
  subset.loci <- named.loci[symbols]
  names(subset.loci) <- NULL
  subset.loci
}


#
# STREAM
#
message ("Began converting eRegulons from nested lists to data.table")
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/stream/"
setwd(work.dir)
getwd()
stream.regionGene <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.qsave"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res <- qs::qread(paste0(work.dir, data.ls[i], "/eRegulons.qsave")) %>% lapply(., "[[", "links") %>%
    lapply(., as.data.frame) %>% Reduce("rbind", .)
  dt <- data.frame(
    Region = paste0(res$seqnames, "-", res$start, "-", res$end), 
    Gene = res$gene, 
    Gene.loc = get_gene_loc(symbols = res$gene, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return (dt[!duplicated(dt) & !is.na(dt$Gene.loc), , drop = FALSE])
})


head(stream.regionGene[[1]])
names(stream.regionGene) <- data.ls
qs::qsave(stream.regionGene, paste0(parent.dir, "stream_region_gene_links.qsave") )
setwd("..")


#
# SCENIC+
#
# Conversion
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scenicplus/"
setwd(work.dir)
scplus.regionGene <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res <- read.csv(paste0(work.dir, data.ls[i], "/eRegulons.csv"), sep = ",")[, c("Region", "Gene")]
  dt <- data.frame(
    Region = res$Region, 
    Gene = res$Gene, 
    Gene.loc = get_gene_loc(symbols = res$Gene, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return (dt[!duplicated(dt) & !is.na(dt$Gene.loc), , drop = FALSE])
})


head(scplus.regionGene)
names(scplus.regionGene) <- data.ls
scplus.regionGene <- scplus.regionGene[!sapply(scplus.regionGene, is.null)]
length(scplus.regionGene)
qs::qsave(scplus.regionGene, paste0(parent.dir, "scenicplus_region_gene_links.qsave") )
setwd("..")


#
# Pando
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/pando/"
setwd(work.dir)
library(Pando)
pando.regionGene <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".qs"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn <- qs::qread(paste0(work.dir, data.ls[i], ".qs")) %>% NetworkModules
  dt <- data.frame(
    Region = grn@meta$regions, 
    Gene = grn@meta$target, 
    Gene.loc = get_gene_loc(symbols = grn@meta$target, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return (dt[!duplicated(dt) & !is.na(dt$Gene.loc), , drop = FALSE])
  message ("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
head(pando.regionGene[[1]])
names(pando.regionGene) <- data.ls
pando.regionGene <- pando.regionGene[!sapply(pando.regionGene, is.null)]
length(pando.regionGene)
qs::qsave(pando.regionGene, paste0(parent.dir, "pando_region_gene_links.qsave"))
setwd("..")


#
# GLUE
#
work.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scglue/"
setwd(work.dir)
scglue.regionGene <- pblapply(seq_along(data.ls), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message ("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn <- read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  gene2peak <- read.csv(paste0(work.dir, data.ls[i], "/gene2peak.csv"), sep = ",")
  colnames(gene2peak) <- c("Gene", "Region", "Weight")
  colnames(grn) <- c("TF", "Gene", "importance")
  res <- gene2peak[gene2peak$Gene %in% grn$Gene,, drop = FALSE]
  dt <- data.frame(
    Region = res$Region, 
    Gene = res$Gene, 
    Gene.loc = get_gene_loc(symbols = res$Gene, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return (dt[!duplicated(dt) | is.na(dt$Gene.loc), , drop = FALSE])
})
head(scglue.regionGene)
names(scglue.regionGene) <- data.ls
scglue.regionGene <- scglue.regionGene[!sapply(scglue.regionGene, is.null)]
length(scglue.regionGene)
qs::qsave(scglue.regionGene, paste0(parent.dir, "scglue_region_gene_links.qsave"))
setwd("..")


# Get list of predicted region-to-gene linkages for each method
method.list <- list(
  stream = stream.regionGene, 
  scenicplus = scplus.regionGene, 
  scglue = scglue.regionGene,
  pando = pando.regionGene
)


################################################################
#                                                              #
#        Benchmark region-gene linkages for each data set      #
#                                                              #
################################################################


# Read correspondence file
top.peaks <- 100
chipseq.hic <- read.csv(file = paste0(parent.dir, "ChIPseq-and-HiC-data.csv"))
dim(chipseq.hic)
head(chipseq.hic)


# Define the function to read Hi-C files
read_hic <- function(fname, org = "hg38", unit = "BP", binsize = 2500000, top.peaks = 100) {
  
  invisible(library(strawr))
  invisible(library(dplyr))
  invisible(ifelse (grepl("^mm", org), chr.ll <- c(1:40, "X", "Y"), chr.ll <- c(1:23, "X", "Y")) )
  messag <- try(straw(norm = "NONE", fname = fname, "1", "1", unit = unit, binsize = binsize))
  if (class(messag) == "try-error") {
    message ("--------> Chromosome names lack a prefix 'chr'")
    chr.ll <- paste0("chr", chr.ll)
  }
  df <- pblapply(chr.ll, function(x) {
    message (x)
    chr.hic <- try(straw(norm = "NONE", fname = fname, x, x, unit = unit, binsize = binsize) %>% arrange(desc(counts)) %>% 
                     head(n = top.peaks))
    if (class(chr.hic) == "try-error") {
      return(data.frame(
        peak1 = NA,
        peak2 = NA,
        counts = NA
      ))
    }
    lapply(1:nrow(chr.hic), function(i) {
      y <- chr.hic[i,]
      data.frame(
        peak1 = paste0(x, "-", format(y$x + 1, scientific = FALSE), "-", format(y$x + binsize, scientific = FALSE) ),
        peak2 = paste0(x, "-", format(y$y + 1, scientific = FALSE), "-", format(y$y + binsize, scientific = FALSE) ), 
        counts = y$counts
      )
    }) %>% Reduce("rbind", .)
  }) %>% na.omit %>% Reduce("rbind", .)
  df <- df[!is.na(df$counts),, drop = FALSE]
  if (!class(messag) == "try-error") {
    df$peak1 <- paste0("chr", df$peak1)
    df$peak2 <- paste0("chr", df$peak2)
  }
  df
}


# Read region-gene pair lists from each Hi-C file
hic.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/Chrom-interaction-validation/"
message ("Begin reading Hi-C file names")
contacts.ll <- pbmclapply(1:nrow(chipseq.hic), mc.cores = detectCores(), function(i) {
  message ("----> Reading Hi-C data for the ", i, "-th dataset: ", chipseq.hic$Dataset[i])
    if (chipseq.hic$Hi.C[i] == "N/A") {
    return(NULL)
  }
  dir.ll <- strsplit(chipseq.hic$Hi.C[i], split = ", ") %>% `[[` (1)
  pblapply(dir.ll, function(x) {
    message("Reading hic files in dir: ", x)
    org <- strsplit(x, split = "-")[[1]][2]
    hic.files <- list.files(paste0(eval.dir, x), pattern = ".hic")
    dir.regionGene <- lapply(hic.files, function(y) {
      read_hic(fname = paste0(eval.dir, x, "/", y), org = org, top.peaks = top.peaks)
    }) %>% Reduce("rbind", .)
  }) %>% Reduce("rbind", .)
})
names(contacts.ll) <- chipseq.hic$Dataset
head(contacts.ll[[1]])
qs::qsave(contacts.ll, paste0(eval.dir, "HiC-contact-matrix-list.qsave"))
message ("Finished reading Hi-C contact matrices")


################################################################
#                                                              #
#                     Region-to-gene calculation               #
#                                                              #
################################################################


# Steps:
# 1. For each dataset, calculate the proportion of region-to-gene linkages against Hi-C contact pairs
# 3. For each dataset, draw box plots for precision, recall, and f scores
# 4. Compare benchmarking methods


# Select data containing more than ten TFs
#tf.covData <- qs::qread(paste0(parent.dir, "data-list-containing-more-than-10-tfs.qsave"))
tf.covData <- chipseq.hic$Dataset
tf.covData


# Define the function to calculate region-to-gene overlaps
# r1 : region-to-gene linkages from regulons saved in data frame
# r2 : bin-to-bin contacts from Hi-C saved in data frame
cal_overlap_regions <- function(r1, r2) {
  
  if (is.null(r1) | is.null(r2)) {
    return(c(precision = 0, recall = 0, f.score = 0))
  }
  r1$Region <- strsplit(r1$Region, split = ";") %>% sapply(., "[", 1)
  r1 <- r1[!is.na(r1$Gene.loc),, drop = FALSE]
  invisible(library(GenomicRanges))
  invisible(library(Signac))
  
  
  forward1 <- findOverlaps(StringToGRanges(r1$Region, sep = c(":", "-")), 
                       StringToGRanges(r2$peak1, sep = c(":", "-")) ) # Overlapping correspondence for <1, 1>
  forward2 <- findOverlaps(StringToGRanges(r1$Gene.loc, sep = c(":", "-")), 
                           StringToGRanges(r2$peak2, sep = c(":", "-")) ) # Overlapping correspondence for <2, 2>
  
  reverse1 <- findOverlaps(StringToGRanges(r1$Region, sep = c(":", "-")), 
                           StringToGRanges(r2$peak2, sep = c(":", "-")) ) # Overlapping correspondence for <1, 2>
  reverse2 <- findOverlaps(StringToGRanges(r1$Gene.loc, sep = c(":", "-")), 
                           StringToGRanges(r2$peak1, sep = c(":", "-")) ) # Overlapping correspondence for <2, 1>
  
  
  invisible(library(tidyverse))
  hits1 <- inner_join(as.data.frame(forward1), as.data.frame(forward2))
  hits2 <- inner_join(as.data.frame(reverse1), as.data.frame(reverse2))
  hits <- rbind(hits1, hits2)
  hits <- hits[!duplicated(hits), , drop = FALSE]
  
  
  prec <- length(unique(hits$queryHits)) / nrow(r1)
  recall <- length(unique(hits$subjectHits)) / nrow(r2)
  c(precision = prec, recall = recall, f.score = 2 * prec * recall / (prec + recall))
}


# Calculate three scores of region-to-gene scores for all methods
library(data.table)
res.ll <- pbmclapply(tf.covData, mc.cores = detectCores(), function(chip.data) {
  message ("----> ChIP-seq data: ", chip.data)
  data.scores <- lapply(seq_along(method.list), function(i) {
    #message (i)
    method <- method.list[[i]]
    if (!chip.data %in% names(method)) {
      return(c(prec = 0, recall = 0, f.score = 0))
    }
    method.scores <- cal_overlap_regions(r1 = method[[chip.data]], 
                                         r2 = contacts.ll[[chip.data]] )
    method.scores[is.na(method.scores)] <- 0.0
    method.scores
  }) %>% Reduce("rbind", .)
  rownames(data.scores) <- names(method.list)
  data.scores
})
length(res.ll)
sapply(res.ll, dim)
names(res.ll) <- tf.covData
regionGene.scores <- res.ll
qs::qsave(res.ll, paste0(parent.dir, "Region-to-gene-benchmark-all.qsave"))
qs::qsave(regionGene.scores, paste0(parent.dir, "Region-to-gene-benchmark.qsave"))


save.image("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/3-benchmark-region-gene.RData")
# load("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/3-benchmark-region-gene.RData")