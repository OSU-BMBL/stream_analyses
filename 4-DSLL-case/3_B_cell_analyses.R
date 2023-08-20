########################################################################
#                                                                      #
# 1. ChIP-seq aggregated signals on B cell type                        #
# 2. Accessibility profiles across cell types and ChIP-seq signals     #
#                                                                      #
########################################################################


# Libraries
library(easypackages)
libs <- c(
  "Seurat",
  "Signac",
  "dplyr", 
  "pbmcapply",
  "parallel"
)
libraries(libs)




# Parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Rdata/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Tables/"
tcf7.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Tables/GSE96199/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
setwd(R.dir)




# source codes
source(paste0(tool.dir, "statistics_tools.R"))
source(paste0(tool.dir, "cistrome_tools.R"))
source(paste0(tool.dir, "transcriptome_tools.R"))
source(paste0(tool.dir, "visual_tools.R"))
source(paste0(tool.dir, "genome_tools.R"))
source(paste0(tool.dir, "IO_tools.R"))




# Load data
cts.en.regs <- qs::qread(paste0(R.dir, "Cell_type_specific_eRegulons.qsave"))
obj <- qs::qread(paste0(R.dir, "lymph_obj_annotated.qsave"))
act.obj <- qs::qread(paste0(R.dir, "Obj_active.qsave"))
dim(obj)
dim(act.obj)
length(cts.en.regs)




# Get box plots of cells, genes, enhancers, and linkages
gene.box <- sapply(cts.en.regs, "[[", "genes") %>% sapply(., length) %>% 
  as.data.frame()
colnames(gene.box) <- "Gene"
gene.box <- cbind(Group = rep(1, nrow(gene.box)), gene.box)

enh.box <- sapply(cts.en.regs, "[[", "enhancers") %>% sapply(., length) %>% 
  as.data.frame()
colnames(enh.box) <- "Enhancer"
enh.box <- cbind(Group = rep(1, nrow(enh.box)), enh.box)

link.box <- sapply(cts.en.regs, "[[", "links") %>% sapply(., length) %>% 
  as.data.frame()
colnames(link.box) <- "Enh-gene"
link.box <- cbind(Group = rep(1, nrow(link.box)), link.box)

p.en.reg.info <- get_boxplot(obj = gene.box, y.lab = "Gene") | 
  get_boxplot(obj = enh.box, y.lab = "Enhancer") | 
  get_boxplot(obj = link.box, y.lab = "Enh-gene")
qs::qsave(p.en.reg.info, paste0(R.dir, "Boxplots_eReghulon_basic_info.qsave"))




# Histogram of # of enhs per genes; th gene per enh
gene.links <- do.call("c", pbmclapply(cts.en.regs, function(x) {
  df <- data.frame(enh = GRangesToString(x$links), 
                   gene = as.data.frame(x$links)[, 6])
  split(df, df$gene) %>% sapply(., nrow)
}, mc.cores = detectCores()))
p.hist.enhs.per.gene <- get_hist(df = gene.links, title = "# of enhs per gene", 
                                 binwidth = 1, ylab = "Frequency")


gene.gr <- CollapseToLongestTranscript(obj[['ATAC']]@annotation)
class(gene.gr)
length(gene.gr)


union.links <- Reduce("rbind", pbmclapply(cts.en.regs, function(x) {
  as.data.frame(x$links)
}, mc.cores = detectCores())) %>% unique
dim(union.links)
union.links <- calc_peak_gene_dist(union.links = union.links, 
                                     gene.gr = gene.gr)
union.links <- cbind(string = GRangesToString(makeGRangesFromDataFrame(df = union.links,
                                                       keep.extra.columns = F)), 
                     union.links)
head(union.links)
rownames(union.links) <- paste0(union.links$string, "_", union.links$gene)
qs::qsave(union.links, paste0(R.dir, "DSLL_enh_gene_distance.qsave"))


enh.ranks <- do.call("c", pblapply(cts.en.regs, function(x) {
  df <- as.data.frame(x$links)
  df <- cbind(string = GRangesToString(makeGRangesFromDataFrame(df = df,
                                                                keep.extra.columns = F)), 
              df)
  rownames(df) <- paste0(df$string, "_", df$gene)
  df <- cbind(df, distance = union.links[rownames(df), "distance"])
  df <- df %>% group_by(gene) %>% arrange(distance, .by_group = T)
  splitted <- split(df, f = df$gene)
  do.call("c", sapply(splitted, function(y) {
    cbind(id = 1:nrow(y), y) %>% pull(id, string)
  }))
}))
p.hist.nth.gene.per.enh <- get_hist(df = enh.ranks, title = "nth. gene per enh", 
                                 binwidth = 1, ylab = "Frequency")
p.hist.gene.enh <- p.hist.enhs.per.gene / p.hist.nth.gene.per.enh
qs::qsave(p.hist.gene.enh, paste0(R.dir, "Histograms_genes_enhs.qsave"))


# B cells
sapply(cts.en.regs, "[[", "cell.type")
normal.en.regs <- cts.en.regs[21:22]
prol.en.regs <- cts.en.regs[39:46]
tumor.en.regs <- cts.en.regs[30:38]
qs::qsave(list(normal.en.regs, prol.en.regs, tumor.en.regs), 
          paste0(R.dir, "B_cell_eReguons.qsave"))
sapply(normal.en.regs, "[[", "TF") 
sapply(prol.en.regs, "[[", "TF") 
sapply(tumor.en.regs, "[[", "TF") 
full.colors[c(9, 8, 6)] 
show_col(full.colors[c(9, 8, 6)])




# Find markers of prol or tumor B cells
Reduce("union", sapply(normal.en.regs, "[[", "de.genes")) 
Reduce("union", sapply(prol.en.regs, "[[", "de.genes")) 
Reduce("union", sapply(tumor.en.regs, "[[", "de.genes"))
setdiff(Reduce("union", sapply(prol.en.regs, "[[", "de.genes")), 
        Reduce("union", sapply(tumor.en.regs, "[[", "de.genes")))
markers <- read.csv(paste0(table.dir, "DSLL_cell_type_markers.csv"))
dim(markers)
markers[, 8:10]




# Tumor prol B cell markers
prol.uniq.genes <- setdiff(Reduce("union", sapply(prol.en.regs, "[[", "de.genes")), 
                           Reduce("union", sapply(tumor.en.regs, "[[", "de.genes")))
prol.uniq.genes
prol.en.grn <- eRegulons_to_eGRN(TFs = lapply(prol.en.regs, "[[", "TF"), 
                                 links = lapply(prol.en.regs, "[[", "links"))
dim(prol.en.grn)
prol.diff.net <- prol.en.grn[prol.en.grn$enhancer %in% 
                             Reduce("union", sapply(prol.en.regs, "[[", "da.enhs")) & 
                             prol.en.grn$gene %in% 
                             Reduce("union", sapply(prol.en.regs, "[[", "de.genes")),]
dim(prol.diff.net)
prol.diff.net[prol.diff.net$gene %in% prol.uniq.genes,]




save.image("3_DSLL_Bcell_analyses.RData")
stop("!!!!!!!")


require(GEOquery)
tcf7.files <- list.files(tcf7.dir, pattern = "_GRCh38.bed")
tcf7.files <- paste0(tcf7.dir, tcf7.files)
tcf7.peaks <- input_peaks(tcf7.files[1])
length(tcf7.peaks) 
head(tcf7.peaks)


tcf7.overlap <- intersect_peaks(x = StringToGRanges(normal.en.regs[[2]]$enhancers), 
                y = tcf7.peaks)
tcf7.overlap <- list(pval = tcf7.overlap$numOverlaps$pval, 
                     overlap = tcf7.overlap$numOverlaps$observed)



# Tumor B cells
setdiff(sapply(tumor.en.regs, "[[", "TF"), 
        sapply(prol.en.regs, "[[", "TF")) 
sapply(prol.en.regs, "[[", "TF") 


# Download ChIP-seq data for EGR1 (PMID: 33980611) (hg38)
# ENCFF637UJN.bed.gz
egr1.peaks <- input_peaks(paste0(table.dir, "ENCFF637UJN.bed"))
egr1.overlap <- intersect_peaks(x = StringToGRanges(tumor.en.regs[[1]]$enhancers), 
                                y = egr1.peaks)
egr1.overlap <- list(pval = egr1.overlap$numOverlaps$pval, 
                     overlap = egr1.overlap$numOverlaps$observed)


# Download ChIP-seq data for STAT3 (hg38)
sapply(tumor.en.regs, "[[", "TF")
stat3.peaks <- input_peaks(paste0(table.dir, "ENCFF029PAA.bed"))
stat3.overlap <- intersect_peaks(x = StringToGRanges(tumor.en.regs[[7]]$enhancers), 
                                y = stat3.peaks)
stat3.overlap <- list(pval = stat3.overlap$numOverlaps$pval, 
                     overlap = stat3.overlap$numOverlaps$observed)
stat3.overlap



# NFIB
sapply(tumor.en.regs, "[[", "TF")
nfic.peaks <- input_peaks(paste0(table.dir, "ENCFF980SDE.bed"))
nfic.overlap <- intersect_peaks(x = StringToGRanges(tumor.en.regs[[4]]$enhancers), 
                                y = nfic.peaks)
nfic.overlap <- list(pval = nfic.overlap$numOverlaps$pval, 
                     overlap = nfic.overlap$numOverlaps$observed)
nfic.overlap

