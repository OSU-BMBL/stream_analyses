
###########################################################
#                                                         #
#                   Global analysis on DSLL               #
#                                                         #
###########################################################


# UMAP plot
# Dotplot-heatmap of TF expression, accessibility, and expression


# Libraries
library(easypackages)
libs <- c("Seurat", "dplyr", "Signac", "pbmcapply", "parallel", "pbapply", 
          "HiCDCPlus", "scales", "ggplot2", "DescTools")
libraries(libs)



# Parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Rdata/"
setwd(R.dir)
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Tables/"



# Functions
source(paste0(tool.dir, "visual_tools.R"))
source(paste0(tool.dir, "transcriptome_tools.R"))
source(paste0(tool.dir, "cistrome_tools.R"))



# Load Rdata
load("6_DSLL_single_gene_analyses.RData")



# Determine TFs, enhancer, and gene
normal.en.grns <- eRegulons_to_eGRN(TFs = sapply(normal.en.regs, "[[", "TF"), 
                                  links = sapply(normal.en.regs, "[[", "links"), 
                                  col = "gene")
tumor.en.grns <- eRegulons_to_eGRN(TFs = sapply(tumor.en.regs, "[[", "TF"), 
                                  links = sapply(tumor.en.regs, "[[", "links"), 
                                  col = "gene")
sapply(normal.en.regs, "[[", "da.enhs")
prol.en.grns
lhx2.enhs <- Reduce("intersect", list(unique(normal.en.grns[normal.en.grns$TF == "LHX2",]$enhancer),
          unique(prol.en.grns[prol.en.grns$TF == "LHX2",]$enhancer),
          unique(tumor.en.grns[tumor.en.grns$TF == "LHX2",]$enhancer)))
length(lhx2.enhs)
range(da.enhs$avg_log2FC)
range(da.enhs$p_val_adj)
lhx2.da.enhs <- da.enhs[da.enhs$cluster == "Tumor B prol" & da.enhs$gene %in% lhx2.enhs & 
          order(da.enhs$avg_log2FC, decreasing = T), "gene"] %>% unique
length(lhx2.da.enhs)
key.enh <- lhx2.da.enhs[5] # marker: ST6GAL1



# Discover the enhancer bound by STAT3, SPI1, and TCF12
cobind.tfs <- c("STAT3", "SPI1", "TCF12")
common.enhs <- Reduce("intersect", list(
  prol.en.grns[prol.en.grns$TF == "STAT3", "enhancer"], 
  prol.en.grns[prol.en.grns$TF == "SPI1", "enhancer"], 
  prol.en.grns[prol.en.grns$TF == "TCF12", "enhancer"]
))
length(common.enhs)
b.da.sig <- b.da.enhs[b.da.enhs$p_val_adj < 0.05 & b.da.enhs$avg_log2FC > 0,]

common.dars <- common.enhs
common.dars.normal <- intersect(common.dars, unique(normal.en.grns$enhancer))
common.dars.prol <- intersect(common.dars, unique(prol.en.grns$enhancer))
common.dars.tumor <- intersect(common.dars, unique(tumor.en.grns$enhancer))
choice.dars <- Reduce("intersect", list(common.dars.normal, common.dars.prol, common.dars.tumor))
prol.degs <- de.genes[de.genes$cluster == "Tumor B prol", "gene"] %>% unique
prol.en.grns[prol.en.grns$enhancer %in% choice.dars & prol.en.grns$gene %in% prol.degs,] # CCDC88C



# Choose key enhancer
key.enh <- "chr14-91366590-91389532"
normal.en.grns[normal.en.grns$enhancer == key.enh, "gene"] # CCDC88C
prol.en.grns[prol.en.grns$enhancer == key.enh, "gene"]
tumor.en.grns[tumor.en.grns$enhancer == key.enh, "gene"]



# TFs
normal.lhx2.tfs <- normal.en.grns[normal.en.grns$enhancer == key.enh, "TF"] %>% unique
prol.lhx2.tfs <- prol.en.grns[prol.en.grns$enhancer == key.enh, "TF"] %>% unique
tumor.lhx2.tfs <- tumor.en.grns[tumor.en.grns$enhancer == key.enh, "TF"] %>% unique




# Select ENCODE datasets
encode.url <- c(
  "https://www.encodeproject.org/files/ENCFF750YNG/@@download/ENCFF750YNG.bigWig", # EGR1
  "https://www.encodeproject.org/files/ENCFF766JTP/@@download/ENCFF766JTP.bigWig", # NFIC
  "https://www.encodeproject.org/files/ENCFF042MFJ/@@download/ENCFF042MFJ.bigWig", # SPI1
  "https://www.encodeproject.org/files/ENCFF011FHF/@@download/ENCFF011FHF.bigWig", # STAT3
  "https://www.encodeproject.org/files/ENCFF579VGX/@@download/ENCFF579VGX.bigWig" # TCF12
)


all.enhs <- Reduce("union", sapply(c(normal.en.regs, prol.en.regs, tumor.en.regs), "[[", "enhancers")) %>% 
  as.data.frame
write.table(all.enhs, paste0(table.dir, "All_enhs.txt"), quote = F, 
            row.names = F, col.names = F)
bed.files <- list(
  "ENCFF492ZRZ.bed", # SPI1
  "ENCFF029PAA.bed", # STAT3
  "ENCFF068YYR.bed" # TCF12
) %>% paste0(table.dir, .) %>% as.list
names(bed.files) <- c("SPI1", "STAT3", "TCF12")


bw.files <- list(
  "../Tables/ENCFF042MFJ.bigWig", # SPI1
  "../Tables/ENCFF011FHF.bigWig", # STAT3
  "../Tables/ENCFF579VGX.bigWig" # TCF12
)

names(bw.files) <- c("SPI1", "STAT3", "TCF12")



# < Here: run Python script >



# Load the profiles
prof.files <- gsub("bigWig", "prof", bw.files)


prof.ls <- pblapply(prof.files, function(x) {
  df <- as.data.frame(read.csv(paste0(table.dir, x), header = F))
  rownames(df) <- all.enhs[, 1]
  df
})
sapply(prof.ls , dim)
names(prof.ls) <- c(
  "SPI1", 
  "STAT3", 
  "TCF12"
)
qs::qsave(prof.ls, paste0(R.dir, "ChIP_seq_profiles.qsave"))
sapply(prol.en.regs, "[[", "TF")
en.regs.enhs <- lapply(prol.en.regs[5:7], "[[", "enhancers")
names(en.regs.enhs) <- names(prof.ls)



# Get profiles of TFs in Tumor B proliferation
prof.tfs <- c(
  "SPI1", 
  "STAT3", 
  "TCF12"
)
reg.chip.profs <- do.call("c", pbmclapply(1:length(prof.tfs), function(i) {
  x <- en.regs.enhs[[prof.tfs[[i]]]] # eRegulon
  ls1 <- lapply(1:i, function(j) {
    if (i != j) {
      x <- intersect(x, en.regs.enhs[[prof.tfs[[j]]]])
    }
    ls2 <- lapply(prof.tfs, function(y) {
      prof.ls[[y]][x,]
    })
    names(ls2) <- prof.tfs
    ls2
  })
  names(ls1) <- paste0(prof.tfs[[i]], "_", prof.tfs[1:i])
  ls1
}, mc.cores = detectCores()))
length(reg.chip.profs)
names(reg.chip.profs)



# Plotting
gg.ls <- c()
tf.colors <- c(
  "#D2001A", 
  "#367E18",
  "#5800FF"
)
names(tf.colors) <- c("TCF12", "SPI1", "STAT3")
for (i in seq_along(reg.chip.profs)) {
  gg.ls[[i]] <- get_prof(df.ls = reg.chip.profs[[i]], norm = "mean",
                         min.max = "local", norm.row = F, 
                         colors = tf.colors, ratio = 1.0)
}
names(gg.ls) <- names(reg.chip.profs)
qs::qsave(gg.ls, paste0(R.dir, "ChIP_seq_profs_on_enhs.qsave"))



# Arrangement of ChIP-seq profiles on enhancers
p.chip.prof <- ggpubr::ggarrange(gg.ls[[1]], ggplot() + theme_void(), ggplot() + theme_void(), 
  gg.ls[[2]], gg.ls[[3]], ggplot() + theme_void(), 
  gg.ls[[4]], gg.ls[[5]], gg.ls[[6]]
)
qs::qsave(p.chip.prof, paste0(R.dir, "Profiles_ChIP_seq_on_enhancers.qsave"))
save_image(p = p.chip.prof, paste0(image.dir, "ChIP_seq_profs_on_enhs.png"), 
            width = 1500, height = 2500)



# Prepare TF binding sites on the key enhancer
load("/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/TFBS_list.rda")
key.tf.sites <- pblapply(prof.tfs, function(x) {
  tf.sites <- TFBS.list$Human$peak[which(TFBS.list$Human$TF == x)]
  subject.query <- findOverlaps(StringToGRanges(key.enh), tf.sites)
  key.sites <- tf.sites[subjectHits(subject.query) %>% unique]
})
names(key.tf.sites) <- prof.tfs
sapply(key.tf.sites, length)
key.gene <- "CCDC88C"
DefaultAssay(b.obj) <- "ATAC"
b.colors <- c("#FFADAD", "#FF5D5D", "#D61C4E")
names(b.colors) <- unique(b.obj$ct.label) %>% rev
show_col(b.colors)
tf.colors <- color.ls[6:8]
tf.colors
colored.tf.sites <- lapply(seq_along(key.tf.sites), function(i) {
  mcols(key.tf.sites[[i]])$color <- rep(tf.colors[[i]], length(key.tf.sites[[i]]))
  key.tf.sites[[i]]
})



# Calculate overlap with Hi-C
# Hi-C: ENCFF053VBX.hic

# ENCFF767OUV
hic.file <- "ENCFF053VBX.hic"
key.enh
construct_features(output_path = paste0(table.dir,"ENCFF053VBX_50kb_GATC"),
                   gen = "Hsapiens", gen_ver = "hg38",
                   sig = "GATC",
                   bin_type = "Bins-uniform",
                   binsize = 50000,
                   chrs = "chr14")


# generate gi_list instance
gi_list <- generate_bintolen_gi_list(
  bintolen_path <- paste0(table.dir,"/ENCFF053VBX_50kb_GATC_bintolen.txt.gz"))


# add .hic counts
gi_list <- add_hic_counts(gi_list, hic_path = paste0(table.dir, "ENCFF767OUV.hic"))


#expand features for modeling
gi_list <- expand_1D_features(gi_list)
set.seed(1010) #HiC-DC downsamples rows for modeling
gi_list <- HiCDCPlus(gi_list) #HiCDCPlus_parallel runs in parallel across ncores
head(gi_list)
qs::qsave(gi_list, paste0(R.dir, "HiC_interact.qsave"))
hic.interact <- as.data.frame(gi_list)
sig.interacts <- hic.interact[hic.interact$chr14.qvalue < 0.05,]
range(sig.interacts$chr14.start1)
hic.interact[1:5, 1:5]



# Get enh-gene linkages
gene.regions <- obj@assays$ATAC@annotation
gene.tss <- resize(get_gene_loc(key.gene = key.gene, gene.regions = gene.regions), 1, 
                   "start")
region.highlight <- sort(resize(do.call("c", colored.tf.sites), 
                               width = 300))
region.highlight2 <- region.highlight[c(2, 5, 6, 7, 8, 13)]
link.gr <- punion(rep(gene.tss, length(region.highlight)),
                  region.highlight, fill.gap = T, ignore.strand = T)
link.gr
mcols(link.gr)$gene <- rep(key.gene, length(region.highlight))
mcols(link.gr)$score <- rep(0.5, length(region.highlight))
mcols(link.gr)$peak <- GRangesToString(region.highlight)
DefaultAssay(b.obj) <- "ATAC"
Links(b.obj) <- link.gr



p.singleGeneCov <- get_coverage_plot(object = b.obj, links = F,
                                     region = paste0("chr14-", min(start(link.gr)), "-", 
                                                     max(end(link.gr))), 
                                     features = key.gene,
                  ranges.group.by = "ct.label", peaks = F, 
                  bigwig = bw.files, region.highlight = NULL, 
                  colors = b.colors, size = 10)
qs::qsave(p.singleGeneCov, paste0(R.dir, "Coverage_plot.qsave"))
save_image(p = p.singleGeneCov, path = paste0(image.dir, "Coverage_plot.png"), 
           width = 3000, height = 800)


# Accessibility profiles of three TFs (Mann-Whitney test p-value)
all.enhs # I still use it
sapply(prol.en.regs, "[[", "TF")
atac.prof <- as.data.frame(read.csv(paste0(table.dir, "lymph_node_lymphoma_14k_atac_cut_sites.prof"), 
                                    header = F))
dim(atac.prof)
rownames(atac.prof) <- all.enhs[, 1]
atac.prof[1:5, 1:5]
reg.atac.profs <- do.call("c", pbmclapply(1:length(prof.tfs), function(i) {
  x <- en.regs.enhs[[prof.tfs[[i]]]] # eRegulon
  ls1 <- lapply(1:i, function(j) {
    if (i != j) {
      x <- intersect(x, en.regs.enhs[[prof.tfs[[j]]]])
    }
    ls2 <- atac.prof[x,]
    ls2
  })
  names(ls1) <- paste0(prof.tfs[[i]], "_", prof.tfs[1:i])
  ls1
}, mc.cores = detectCores()))
length(reg.atac.profs)
names(reg.atac.profs)
dim(reg.atac.profs[[1]])



# Plot the profile
require(DescTools)
names(reg.atac.profs) <- c(
  "SPI1", 
  "SPI1+STAT3", 
  "STAT3", 
  "SPI1+TCF12", 
  "STAT3+TCF12", 
  "TCF12")
names(reg.atac.profs)
tf.colors
comb.colors <- c(tf.colors[[2]], MixColor(tf.colors[[2]], 
                 tf.colors[[3]], amount1 = 0.5), tf.colors[[3]], 
                 MixColor(tf.colors[[1]], 
                          tf.colors[[2]], amount1 = 0.5), 
                 MixColor(tf.colors[[1]], 
                          tf.colors[[3]], amount1 = 0.5), tf.colors[[1]])
names(comb.colors) <- names(reg.atac.profs)
show_col(comb.colors)
p.atac.prof <- get_prof(df.ls = reg.atac.profs, colors = comb.colors,
                        norm = "mean", ratio = 1.0,
                        min.max = "global", norm.row = F, 
                        y.lab = "Accessibility in Tumor B prol",
                        legend.position = "top")
qs::qsave(p.atac.prof, paste0(R.dir, "ATAC_prof_key_enh.qsave"))
save_image(p = p.atac.prof, path = paste0(image.dir, "ATAC_prof_key_enh.png"), 
           width = 1500, height = 1000)


# Single gene (not SPI1, STAT3, and TCF12)
single.tfs <- names(color.ls)[2:8]
single.tf.sites <- pblapply(single.tfs, function(x) {
  tf.sites <- TFBS.list$Human$peak[which(TFBS.list$Human$TF == x)]
  subject.query <- findOverlaps(StringToGRanges(key.enh), tf.sites)
  key.sites <- tf.sites[subjectHits(subject.query) %>% unique]
})
names(single.tf.sites) <- single.tfs
sapply(single.tf.sites, length)
single.gene <- "CCDC88C"
single.colors <- color.ls[2:8]
single.colors
colored.single <- lapply(seq_along(single.tf.sites), function(i) {
  mcols(single.tf.sites[[i]])$color <- rep(single.colors[[i]], 
                                           length(single.tf.sites[[i]]))
  single.tf.sites[[i]]
})
head(colored.single)
tail(colored.single)



# Get enh-gene linkages for single gene
single.region.highlight <- sort(resize(do.call("c", colored.single), 
                                width = 300))
single.link.gr <- punion(rep(gene.tss, length(single.region.highlight)),
                  single.region.highlight, fill.gap = T, ignore.strand = T)
single.link.gr
mcols(single.link.gr)$gene <- rep(key.gene, length(single.region.highlight))
mcols(single.link.gr)$score <- rep(0.5, length(single.region.highlight))
mcols(single.link.gr)$peak <- GRangesToString(single.region.highlight)
DefaultAssay(b.obj) <- "ATAC"
Links(b.obj) <- single.link.gr
single.region.highlight2 <- single.region.highlight[c(1, 2, 5, 8, 10, 17, 22)]



p.singleGeneMultiTFsCov <- get_coverage_plot(object = b.obj, links = F,
                                     region = paste0("chr14-", min(start(single.link.gr)), "-", 
                                                     max(end(single.link.gr))), 
                                     features = key.gene,
                                     ranges.group.by = "ct.label", peaks = T, 
                                     region.highlight = single.region.highlight2, 
                                     colors = b.colors, size = 10)
qs::qsave(p.singleGeneMultiTFsCov, paste0(R.dir, "single_gene_multi_TFs_coverage_plot.qsave"))
save_image(p = p.singleGeneMultiTFsCov, 
           path = paste0(image.dir, "single_gene_multi_TFs_coverage_plot.png"), 
           width = 3000, height = 1500)



save.image("6_DSLL_single_gene_analyses.RData")
