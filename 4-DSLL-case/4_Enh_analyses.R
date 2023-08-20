#############################################################
#                                                           #
# 1. Inspect selected TFs from three cell types             #
# 2. For constant TF, check the changes of relationships    #
#    and enhancers                                          #
# 3. For unique TFs, show the core enhancers                #
# 4. Visualize the results using networks, coverage, track, #
#    and pathways                                           #
# 5. ChIP-seq and Hi-C validation                           #
#                                                           #
#############################################################



# Parameters
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Images/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"



# Sources
source(paste0(tool.dir, "cistrome_tools.R"))
source(paste0(tool.dir, "IO_tools.R"))
source(paste0(tool.dir, "genome_tools.R"))




# Construct eGRNs using eRegulons for normal B cells
normal.uniq.genes <- Reduce("union", sapply(normal.en.regs, "[[", "de.genes"))
normal.en.grn <- eRegulons_to_eGRN(TFs = lapply(normal.en.regs, "[[", "TF"), 
                                 links = lapply(normal.en.regs, "[[", "links"))
dim(normal.en.grn)
normal.diff.net <- normal.en.grn[normal.en.grn$enhancer %in% 
                               Reduce("union", sapply(normal.en.regs, "[[", "da.enhs")) & 
                               normal.en.grn$gene %in% 
                               Reduce("union", sapply(normal.en.regs, "[[", "de.genes")),]
dim(normal.diff.net)
normal.diff.net[normal.diff.net$gene %in% normal.uniq.genes,]




# Construct eGRNs using eRegulons for prol B cells
prol.uniq.genes <- Reduce("union", sapply(prol.en.regs, "[[", "de.genes"))
prol.en.grn <- eRegulons_to_eGRN(TFs = lapply(prol.en.regs, "[[", "TF"), 
                                   links = lapply(prol.en.regs, "[[", "links"))
dim(prol.en.grn)
prol.diff.net <- prol.en.grn[prol.en.grn$enhancer %in% 
                                   Reduce("union", sapply(prol.en.regs, "[[", "da.enhs")) & 
                                   prol.en.grn$gene %in% 
                                   Reduce("union", sapply(prol.en.regs, "[[", "de.genes")),]
dim(prol.diff.net)
prol.diff.net[prol.diff.net$gene %in% prol.uniq.genes,]




# Construct eGRNs using eRegulons for tumor B cells
tumor.uniq.genes <- Reduce("union", sapply(tumor.en.regs, "[[", "de.genes"))
tumor.en.grn <- eRegulons_to_eGRN(TFs = lapply(tumor.en.regs, "[[", "TF"), 
                                   links = lapply(tumor.en.regs, "[[", "links"))
dim(tumor.en.grn)
tumor.diff.net <- tumor.en.grn[tumor.en.grn$enhancer %in% 
                                   Reduce("union", sapply(tumor.en.regs, "[[", "da.enhs")) & 
                                   tumor.en.grn$gene %in% 
                                   Reduce("union", sapply(tumor.en.regs, "[[", "de.genes")),]
dim(tumor.diff.net)
tumor.diff.net[tumor.diff.net$gene %in% tumor.uniq.genes,]



# Constant TF LHX2 across three cell types
normal.lhx2 <- normal.diff.net[normal.diff.net$TF == "LHX2",]
prol.lhx2 <- prol.diff.net[prol.diff.net$TF == "LHX2",]
tumor.lhx2 <- tumor.diff.net[tumor.diff.net$TF == "LHX2",]




# Co-factors of LHX2
normal.cof.lhx2 <- normal.diff.net[normal.diff.net$enhancer %in% 
                                     normal.lhx2$enhancer,]
prol.cof.lhx2 <- prol.diff.net[prol.diff.net$enhancer %in% 
                                 prol.cof.lhx2$enhancer,]
tumor.cof.lhx2 <- tumor.diff.net[tumor.diff.net$enhancer %in% 
                                   tumor.lhx2$enhancer,]
require(RCy3)
require(igraph)
core.tfs <- Reduce("union", list(unique(normal.diff.net[normal.diff.net$gene %in% normal.uniq.genes,]$TF), 
                     unique(prol.diff.net[prol.diff.net$gene %in% prol.uniq.genes,]$TF), 
                     unique(tumor.diff.net[tumor.diff.net$gene %in% tumor.uniq.genes,]$TF)))
length(core.tfs)
tf.colors <- c(
  "#FC3C3C", 
  "#5463FF", 
  "#32E0C4", 
  "#FF5B00", 
  "#9818D6",
  "#C40B13",
  "#EA2BA2", 
  "#0F00FF", 
  "#6AF79A",
  "#FFDD00"
)
names(tf.colors) <- core.tfs
show_col(tf.colors)


# LHX2 eGRN in normal B cells
normal.cof.lhx2
colnames(normal.cof.lhx2)
df_to_igraph_net(df = normal.cof.lhx2, color.df = tf.colors, 
                 path = paste0(image.dir, "Normal_B_LHX2.png"))


# LHX2 eGRN in prol B cells
prol.cof.lhx2
colnames(prol.cof.lhx2)
df_to_igraph_net(df = prol.cof.lhx2, color.df = tf.colors, 
                 path = paste0(image.dir, "Prol_B_LHX2.png"))


# LHX2 enhancer accessibility
b.obj <- subset(act.obj, subset = ct.label == c("Normal B", "Tumor B prol", "Tumor B"))
dim(b.obj)
unique(normal.cof.lhx2$enhancer)
unique(prol.cof.lhx2$enhancer)
unique(tumor.cof.lhx2$enhancer)
DefaultAssay(b.obj) <- "ATAC"
prol.enhs <- setdiff(unique(prol.cof.lhx2$enhancer), 
                     unique(normal.cof.lhx2$enhancer))
prol.markers <- setdiff(unique(prol.cof.lhx2$gene), 
                        unique(normal.cof.lhx2$gene))[-2]


prol.core.enhs <- unique(prol.cof.lhx2[prol.cof.lhx2$enhancer %in% prol.enhs & 
                       prol.cof.lhx2$gene %in% prol.markers, "enhancer"])
Idents(b.obj) <- b.obj$ct.label
p.vln.prol.lhx2 <- VlnPlot(b.obj, features = prol.core.enhs, group.by = "ct.label", 
                           fill.by = "ct.label", ncol = 4, 
        idents = c("Normal B", "Tumor B prol"), cols = c("#357C3C", "#FF0075"))
qs::qsave(p.vln.prol.lhx2, paste0(R.dir, "VlnPlot_prol_LHX2.qsave"))


prol.lhx2.tfs <- unique(prol.cof.lhx2[prol.cof.lhx2$enhancer %in% prol.enhs & 
                                 prol.cof.lhx2$gene %in% prol.markers, "TF"])




# NFIC data from ENCODE (ENCFF980SDE.bed)
require(Signac)
sapply(prol.en.regs, "[[", "TF")
nfic.peaks <- input_peaks(paste0(table.dir, "ENCFF980SDE.bed"))
spi1.peaks <- input_peaks(paste0(table.dir, "ENCFF492ZRZ.bed"))
stat3.peaks <- input_peaks(paste0(table.dir, "ENCFF029PAA.bed"))
tcf12.peaks <- input_peaks(paste0(table.dir, "ENCFF897RYA.bed"))
prol.chipseq <- list(
  nfic = nfic.peaks,
  spi1 = spi1.peaks,
  stat3 = stat3.peaks,
  tcf12 = tcf12.peaks
)
qs::qsave(prol.chipseq, paste0(R.dir, "Prol_ChIP_seq_peaks.qsave"))


# Calculate overlap
require(pbapply)
require(pbmcapply)
prol.overlap.m <- rbindlist(pbmclapply(prol.en.regs[4:7], function(x) {
  lapply(prol.chipseq, function(y) {
    overlap <- intersect_peaks(x = StringToGRanges(x$enhancers), 
                    y = y, alternative = "greater", n.times = 100)
    return(overlap$numOverlaps$zscore)
  })
}, mc.cores = 4))
prol.overlap.m <- as.matrix(prol.overlap.m)
rownames(prol.overlap.m) <- sapply(prol.en.regs, "[[", "TF")[4:7]
colnames(prol.overlap.m) <- sapply(prol.en.regs, "[[", "TF")[4:7]
qs::qsave(prol.overlap.m, paste0(R.dir, "Prol_overlap_matrix.qsave"))
require(ComplexHeatmap)
require(circlize)
range(prol.overlap.m)
col_fun = colorRamp2(c(range(prol.overlap.m)[1], 
                       round(sum(range(prol.overlap.m)) / 2), 
                       range(prol.overlap.m)[2]), c("#21209C", "#FF5F00", "yellow"))
p.heatmap.prol.chipseq <- Heatmap(as.matrix(prol.overlap.m), show_row_dend = F,
        row_title = NULL, name = "Overlap z-score", cluster_rows = F,
        show_column_dend = F, cluster_columns = F, col = col_fun, 
        rect_gp = gpar(col = NA, lwd = 0.5))
qs::qsave(p.heatmap.prol.chipseq, paste0(R.dir, "Heatmap_overlap_ChIP_seq.qsave"))
# Hi-C data was used for benchmarking only in SCENIC+




# egr1 eGRN in tumor B cells
tumor.egr1 <- tumor.diff.net[tumor.diff.net$TF == "EGR1",]
tumor.cof.egr1 <- tumor.diff.net[tumor.diff.net$enhancer %in% 
                                   tumor.egr1$enhancer,]
tumor.cof.egr1
colnames(tumor.cof.egr1)
df_to_igraph_net(df = tumor.cof.egr1, color.df = tf.colors, 
                 path = paste0(image.dir, "Tumor_B_EGR1.png"))
unique(tumor.cof.egr1$TF)
prol.cof.egr1 <- prol.diff.net[prol.diff.net$TF %in% unique(tumor.cof.egr1$TF) & 
                                 prol.diff.net$enhancer %in% 
                                 unique(tumor.cof.egr1$enhancer),]





tumor.enhs <- setdiff(unique(tumor.cof.egr1$enhancer), 
                     unique(prol.cof.egr1$enhancer))
tumor.markers <- setdiff(unique(tumor.cof.egr1$gene), 
                        unique(prol.cof.egr1$gene))
tumor.core.enhs <- unique(tumor.cof.egr1[tumor.cof.egr1$enhancer %in% tumor.enhs & 
                                         tumor.cof.egr1$gene %in% tumor.markers, "enhancer"])
p.vln.tumor.egr1 <- VlnPlot(b.obj, features = tumor.core.enhs, group.by = "ct.label", 
                           fill.by = "ct.label", 
                           idents = c("Tumor B prol", "Tumor B"), cols = c("#187498", "#357C3C"))
qs::qsave(p.vln.tumor.egr1, paste0(R.dir, "VlnPlot_tumor_EGR1.qsave"))
df_to_igraph_net(df = prol.cof.egr1, color.df = tf.colors, 
                 path = paste0(image.dir, "Prol_B_EGR1.png"))




tumor.egr1.tfs <- unique(tumor.cof.egr1[tumor.cof.egr1$enhancer %in% tumor.enhs & 
                                        tumor.cof.egr1$gene %in% tumor.markers, "TF"])


egr1.peaks <- input_peaks(paste0(table.dir, "ENCFF637UJN.bed"))
tumor.chipseq <- c(egr1 = egr1.peaks, prol.chipseq)
names(tumor.chipseq)




# Calculate overlap
sapply(tumor.en.regs[c(1, 5:8)], "[[", "TF")
tumor.overlap.m <- rbindlist(pbmclapply(tumor.en.regs[c(1, 5:8)], function(x) {
  lapply(tumor.chipseq, function(y) {
    overlap <- intersect_peaks(x = StringToGRanges(x$enhancers), 
                               y = y, alternative = "greater", n.times = 100)
    return(overlap$numOverlaps$zscore)
  })
}, mc.cores = 4))
tumor.overlap.m <- as.matrix(tumor.overlap.m)
rownames(tumor.overlap.m) <- sapply(tumor.en.regs, "[[", "TF")[c(1, 5:8)]
colnames(tumor.overlap.m) <- sapply(tumor.en.regs, "[[", "TF")[c(1, 5:8)]
qs::qsave(tumor.overlap.m, paste0(R.dir, "tumor_overlap_matrix.qsave"))
require(ComplexHeatmap)
require(circlize)
range(tumor.overlap.m)
col_fun = colorRamp2(c(range(tumor.overlap.m)[1], 
                       round(sum(range(tumor.overlap.m)) / 2), 
                       range(tumor.overlap.m)[2]), c("#21209C", "#FF5F00", "yellow"))
p.heatmap.tumor.chipseq <- Heatmap(as.matrix(tumor.overlap.m), show_row_dend = F,
                                  row_title = NULL, name = "Overlap z-score", cluster_rows = F,
                                  show_column_dend = F, cluster_columns = F, col = col_fun, 
                                  rect_gp = gpar(col = NA, lwd = 0.5))
qs::qsave(p.heatmap.tumor.chipseq, paste0(R.dir, "Heatmap_overlap_ChIP_seq.qsave"))




save.image("4_DSLL_enh_analyses.RData")
