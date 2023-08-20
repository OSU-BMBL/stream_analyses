###########################################################
#                                                         #
#                   Global analysis on DSLL               #
#                                                         #
###########################################################


# Libraries
library(easypackages)
libs <- c("Seurat", "dplyr", "Signac", "pbmcapply", "parallel")
libraries(libs)



# Parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Rdata/"
setwd(R.dir)
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Tables/"
scratch.dir <- "/fs/ess/scratch/PCON0022/liyang/stream/Case_1_DSLL_cell_type/backup_rda/"


load("5_DSLL_global_analyses.RData")



# Load data
act.obj <- qs::qread(paste0(R.dir, "Obj_active.qsave"))
dim(act.obj)
levels(act.obj$ct.label)
act.obj$ct.label <- factor(act.obj$ct.label, 
                           levels = levels(act.obj$ct.label)[c(5, 7, 1, 10, 4, 2, 
                                                               9, 8, 6, 3, 11)])
write.csv(table(act.obj$ct.label), paste0(table.dir, "Cell_type_cells.csv"), 
          quote = F)
levels(act.obj$ct.label)
col.list <- c("#5800FF", "#0096FF", "#72FFFF", 
              "#A63EC5", "#CE49BF", "#ECC5FB", 
              "#FFADAD", "#FF5D5D", "#D61C4E", 
              "#EABF9F", "#F0E9D2")
show_col(col.list)
p.act.umap <- get_simple_UMAP(obj = act.obj, reduction = "umap.rna", 
                              group.by = "ct.label", 
                              cols = setNames(col.list, 
                                levels(act.obj$ct.label)))
p.act.umap
qs::qsave(p.act.umap, paste0(R.dir, "UMAP_celltype_act.qsave"))
save_image(p.act.umap, paste0(image.dir, "UMAP_active_cells.png"), 
           width = 1500, height = 1500)



# Dotplot-heatmap of TF expression, accessibility, expression, 
# or cell type enrichment p-values
cts.en.regs <- qs::qread(paste0(R.dir, "Cell_type_specific_eRegulons.qsave"))
length(cts.en.regs)
sapply(cts.en.regs, "[[", "cell.type")
cts.en.regs <- cts.en.regs[c(
  7:14, 
  1:6, 15:20, 
  21:22, 39:46, 30:38
)]
sapply(cts.en.regs, "[[", "cell.type")
cts.tfs <- sapply(cts.en.regs, "[[", "TF")
AverageExpression(act.obj, assay = "RNA", features = cts.tfs, group.by = "ct.label")$RNA




# DEGs and DARs
de.genes <- qs::qread(paste0(R.dir, "DEGs.qsave"))
da.enhs <- qs::qread(paste0(R.dir, "DARs.qsave"))
dim(de.genes)
dim(da.enhs)
de.genes <- de.genes[de.genes$p_val_adj < 0.05 & de.genes$avg_log2FC > 0,]
da.enhs <- da.enhs[da.enhs$p_val_adj < 0.05 & da.enhs$avg_log2FC > 0,]
dim(de.genes)
dim(da.enhs)
cts.en.regs <- pbmclapply(cts.en.regs, function(x) {
  gene.df <- de.genes[de.genes$cluster == x$cell.type & de.genes$gene %in% x$genes,]
  enh.df <- da.enhs[da.enhs$cluster == x$cell.type & da.enhs$gene %in% x$enhancers,]
  x$de.genes <- head(gene.df$gene, n = Inf)
  x$da.enhs <- head(enh.df$gene, n = Inf)
  if (x$TF %in% de.genes[de.genes$cluster == x$cell.type & de.genes$p_val_adj < 0.05 & 
                         de.genes$avg_log2FC > 0,]$gene) {
    x$de.TF <- T
  } else {
    x$de.TF <- F
  }
  x
}, mc.cores = detectCores())
sapply(cts.en.regs, "[[", "de.genes") %>% sapply(., length)
sapply(cts.en.regs, "[[", "da.enhs") %>% sapply(., length)



# Heatmap of DEGs
act.obj <- ScaleData(act.obj)
exp.m <- act.obj@assays$RNA@data
top10.exp <- Reduce("rbind", pbmclapply(cts.en.regs, function(x) {
  AverageExpression(act.obj, assays = "RNA", features = x$de.genes, 
                    group.by = "ct.label")$RNA %>% apply(., 1, scale) %>% 
    t %>% apply(., 2, mean)
}))
rownames(top10.exp) <- seq_along(cts.en.regs)
colnames(top10.exp) <- levels(act.obj$ct.label)
heatmap(t(scale(t(top10.exp))), scale = "row", Rowv = NA, Colv = NA)



# Expression including all genes
all.exp <- Reduce("rbind", pbmclapply(cts.en.regs, function(x) {
  AverageExpression(act.obj, assays = "RNA", features = x$genes, 
                    group.by = "ct.label")$RNA %>% apply(., 1, scale) %>% 
    t %>% apply(., 2, mean)
}))
rownames(all.exp) <- seq_along(cts.en.regs)
colnames(all.exp) <- levels(act.obj$ct.label)
heatmap(t(scale(t(all.exp))), scale = "row", Rowv = NA, Colv = NA)



# Heatmap of DARs
DefaultAssay(act.obj) <- "ATAC"
act.obj <- ScaleData(act.obj)
acc.m <- act.obj@assays$ATAC@data
top10.acc <- Reduce("rbind", pblapply(cts.en.regs[which(sapply(cts.en.regs, 
                                                               "[[", "da.enhs") %>% 
                                                          sapply(., length) >= 0)], 
                                      function(x) {
                                        if (length(x$da.enhs) > 0) {
                                          AverageExpression(act.obj, assays = "ATAC", features = x$da.enhs, 
                                                            group.by = "ct.label")$ATAC %>% apply(., 1, scale) %>% 
                                            t %>% apply(., 2, mean)
                                        } else {
                                          rep(0, length(unique(act.obj$ct.label)))
                                        }
}))
rownames(top10.acc) <- which(sapply(cts.en.regs, "[[", "da.enhs") %>% 
                               sapply(., length) >= 0)
colnames(top10.acc) <- levels(act.obj$ct.label)
heatmap(t(scale(t(top10.acc))), scale = "row", Rowv = NA, Colv = NA)



# ATAC including all peaks
all.acc <- Reduce("rbind", pblapply(cts.en.regs[which(sapply(cts.en.regs, 
                                                               "[[", "enhancers") %>% 
                                                          sapply(., length) > 0)], 
                                      function(x) {
                                        AverageExpression(act.obj, assays = "ATAC", features = x$enhancers, 
                                                          group.by = "ct.label")$ATAC %>% apply(., 1, scale) %>% 
                                          t %>% apply(., 2, mean)
                                      }))
rownames(all.acc) <- which(sapply(cts.en.regs, "[[", "enhancers") %>% 
                               sapply(., length) > 0)
colnames(all.acc) <- levels(act.obj$ct.label)
heatmap(t(scale(t(all.acc))), scale = "row", Rowv = NA, Colv = NA)
# ATAC is not very obvious but acceptable; I will showcase it via dotplot


range(top10.acc)
require(matrixStats)
acc.rank <- rowRanks(top10.acc)
rownames(acc.rank) <- rownames(top10.acc)
colnames(acc.rank) <- colnames(top10.acc)



range(all.acc)
require(matrixStats)
all.acc.rank <- rowRanks(all.acc)
rownames(all.acc.rank) <- rownames(all.acc)
colnames(all.acc.rank) <- colnames(all.acc)



# Prepare to generate real dotplot-heatmap using top10.exp and acc.rank
en.reg.ids <- which(sapply(cts.en.regs, 
                           "[[", "de.genes") %>% 
                      sapply(., length) > 0)
length(en.reg.ids)
dim(top10.exp)
dim(acc.rank)
core.en.regs <- cts.en.regs[en.reg.ids]
core.exp <- top10.exp[en.reg.ids,]
core.acc <- apply(acc.rank, 1, max) + 1 - acc.rank
core.acc[apply(core.acc, 1, sd) == 0,] <- max(core.acc)
dim(core.exp)
dim(core.acc)
qs::qsave(core.en.regs, paste0(R.dir, "Core_eRegulons.qsave"))



# Generate dotplot-heatmap
sapply(core.en.regs, "[[", "cell.type")
exp.acc.summ <- as(core.exp, "sparseMatrix") %>% summary
colnames(exp.acc.summ) <- c("rowv", "columnv", "EX z-score")
head(exp.acc.summ)
new.col <- sapply(1:nrow(exp.acc.summ), function(i) {
  core.acc[exp.acc.summ$rowv[i], exp.acc.summ$columnv[i]]
})
exp.acc.summ <- cbind(exp.acc.summ, "AC rank" = new.col)
exp.acc.summ <- cbind(exp.acc.summ, label = cut(exp.acc.summ$`AC rank`, 
                                               breaks = 5))
levels(exp.acc.summ$label) <- seq(3, 11, by = 2)
exp.acc.summ <- cbind(exp.acc.summ, size = exp.acc.summ$label)
levels(exp.acc.summ$size) <- 5:1
exp.acc.summ[1:5,]
red.grad <- colorRampPalette(c("#FFAE6D", "#D2001A"))(8)
orange.grad <- colorRampPalette(c("#FFD8A9", "#FFAE6D"))(4)
mid.col <- "#F8F9D7"
blue.grad <- colorRampPalette(c("#0096FF", "#E8F9FD"))(2)
ht.colours <- c(blue.grad, mid.col, orange.grad, red.grad) %>% unique
show_col(ht.colours)
ht.breaks <- c(-0.6, 2.4)
row.labs <- paste0("eR", sapply(seq_along(core.en.regs), function(i) {
  paste0(i, "_", core.en.regs[[i]]$TF)
}))
row.anno <- sapply(core.en.regs, function(x) {
  x$cell.type
})
exp.acc.summ$rowv <- factor(row.labs[exp.acc.summ$rowv], levels = row.labs)
exp.acc.summ$columnv <- factor(levels(act.obj$ct.label)[exp.acc.summ$columnv], 
                               levels = levels(act.obj$ct.label))
p.dot.heat.ex.ac <- get_dotplot_heatmap(df = exp.acc.summ, ht.colours = ht.colours, 
                                        ht.fill = "EX z-score", ht.breaks = ht.breaks, 
                                        dot.name = "AC rank", dot.size = "size", 
                                        dot.breaks = "label", 
                                        row.anno = setNames(col.list, 
                                                            levels(act.obj$ct.label))[sapply(core.en.regs, 
                                                                                             "[[", "cell.type")], 
                                        hjust = 1, vjust = 1, angle = 30)
qs::qsave(p.dot.heat.ex.ac, paste0(R.dir, "Dotplot_heatmap_EX_AC.qsave"))
save_image(p.dot.heat.ex.ac, paste0(image.dir, "Dotplot_heatmap_EX_AC.png"), 
           width = 2400, height = 3600)



# Generate dotplot-heatmap for all genes and peaks
sapply(cts.en.regs, "[[", "cell.type")
all.exp.acc.summ <- as(all.exp, "sparseMatrix") %>% summary
colnames(all.exp.acc.summ) <- c("rowv", "columnv", "EX z-score")
head(all.exp.acc.summ)
all.new.col <- sapply(1:nrow(all.exp.acc.summ), function(i) {
  all.acc.rank[all.exp.acc.summ$rowv[i], all.exp.acc.summ$columnv[i]]
})
all.exp.acc.summ <- cbind(all.exp.acc.summ, "AC rank" = all.new.col)
all.exp.acc.summ <- cbind(all.exp.acc.summ, label = cut(all.exp.acc.summ$`AC rank`, 
                                                breaks = 5))
levels(all.exp.acc.summ$label) <- seq(3, 11, by = 2)
all.exp.acc.summ <- cbind(all.exp.acc.summ, size = all.exp.acc.summ$label)
levels(all.exp.acc.summ$size) <- 5:1
all.exp.acc.summ[1:5,]
red.grad <- colorRampPalette(c("#FFAE6D", "#D2001A"))(8)
orange.grad <- colorRampPalette(c("#FFD8A9", "#FFAE6D"))(4)
mid.col <- "#F8F9D7"
blue.grad <- colorRampPalette(c("#0096FF", "#E8F9FD"))(2)
ht.colours <- c(blue.grad, mid.col, orange.grad, red.grad) %>% unique
show_col(ht.colours)
all.ht.breaks <- c(-0.3, 1.7)
all.row.labs <- paste0("eR", sapply(seq_along(cts.en.regs), function(i) {
  paste0(i, "_", cts.en.regs[[i]]$TF)
}))
all.exp.acc.summ$rowv <- factor(all.row.labs[all.exp.acc.summ$rowv], levels = all.row.labs)
all.exp.acc.summ$columnv <- factor(levels(act.obj$ct.label)[all.exp.acc.summ$columnv], 
                               levels = levels(act.obj$ct.label))
p.dot.heat.ex.ac.all <- get_dotplot_heatmap(df = all.exp.acc.summ, ht.colours = ht.colours, 
                                        ht.fill = "EX z-score", ht.breaks = all.ht.breaks, 
                                        dot.name = "AC rank", dot.size = "size", 
                                        dot.breaks = "label", 
                                        row.anno = setNames(col.list, 
                                                            levels(act.obj$ct.label))[sapply(cts.en.regs, 
                                                                                             "[[", "cell.type")],
                                        hjust = 1, vjust = 1, angle = 30)
qs::qsave(p.dot.heat.ex.ac.all, paste0(R.dir, "Dotplot_heatmap_EX_AC_all.qsave"))
save_image(p.dot.heat.ex.ac.all, paste0(image.dir, "Dotplot_heatmap_EX_AC_all.png"), 
           width = 2000, height = 3200)



# Histogram of # of enhs per genes; th gene per enh
gene.links <- do.call("c", pbmclapply(core.en.regs, function(x) {
  df <- data.frame(enh = GRangesToString(x$links), 
                   gene = as.data.frame(x$links)[, 6])
  split(df, df$gene) %>% sapply(., nrow)
}, mc.cores = detectCores()))



gene.gr <- CollapseToLongestTranscript(obj[['ATAC']]@annotation)
class(gene.gr)
length(gene.gr)


union.links <- Reduce("rbind", pbmclapply(core.en.regs, function(x) {
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


enh.ranks <- do.call("c", pblapply(seq_along(core.en.regs), function(i) {
  message (i)
  x <- core.en.regs[[i]]
  df <- as.data.frame(x$links)
  df <- cbind(string = GRangesToString(makeGRangesFromDataFrame(df = df,
                                                                keep.extra.columns = F)), 
              df)
  rownames(df) <- paste0(df$string, "_", df$gene)
  df <- cbind(df, distance = union.links[rownames(df), "distance"])
  df <- df %>% group_by(gene) %>% dplyr::arrange(distance, .by_group = T)
  splitted <- split(df, f = df$gene)
  xx <- sapply(splitted, function(y) {
    cbind(id = 1:nrow(y), y) %>% pull(id, string)
  })
  
  if (is.list(xx)) {
    return(do.call("c", xx))
  } else {
    return(xx)
  }
}))

hist.df <- rbind(cbind(rep("# of enhs per gene", length(gene.links)),
                       gene.links), 
                 cbind(rep("nth. gene per enh.", length(enh.ranks)),
                       enh.ranks)) %>% as.data.frame
hist.df[, 2] <- as.numeric(hist.df[, 2])
p.hist.gene.enh <- get_multi_hist(df = hist.df, ylab = "Frequency")
qs::qsave(p.hist.gene.enh, paste0(R.dir, "Histograms_genes_enhs.qsave"))
save_image(p.hist.gene.enh, paste0(image.dir, "Histograms_gene_links_enh_ranks.png"), 
           width = 2000, height = 1500)



# Get fraction of overlaps
length(core.en.regs)
sapply(core.en.regs, "[[", "cell.type")
core.dar.overlap <- Reduce("rbind", pbmclapply(core.en.regs, function(x) {
  xx <- length(x$da.enhs)
  sapply(core.en.regs, function(y) {
    yy <- length(y$da.enhs)
    overlap <- length(intersect(x$da.enhs, y$da.enhs))
    overlap / xx
  })
}, mc.cores = detectCores()))
dim(core.dar.overlap)
rownames(core.dar.overlap) <- paste0("eR", seq_along(core.en.regs), "_", 
                                     sapply(core.en.regs, "[[", "TF"))
colnames(core.dar.overlap) <- paste0("eR", seq_along(core.en.regs), "_", 
                                     sapply(core.en.regs, "[[", "TF"))
core.dar.overlap[1:3, 1:3]

core.enh.overlap <- Reduce("rbind", pbmclapply(core.en.regs, function(x) {
  xx <- length(x$enhancers)
  sapply(core.en.regs, function(y) {
    yy <- length(y$enhancers)
    overlap <- length(intersect(x$enhancers, y$enhancers))
    overlap / xx
  })
}, mc.cores = detectCores()))
dim(core.enh.overlap)
rownames(core.enh.overlap) <- paste0("eR", seq_along(core.en.regs), "_", 
                                     sapply(core.en.regs, "[[", "TF"))
core.enh.overlap[1:3, 1:3]
p.heatmap.fr.all <- get_simple_heatmap(mat = core.enh.overlap, splitted = NULL, name = "fr.overlap", 
                   bottom = 0.0, interium = 0.5, up = 1.0, labels = c("0.0", "1.0"),
                   at = c(0, 1), title_position = "leftcenter-rot", 
                   row.name.colors = setNames(col.list, 
                                              levels(act.obj$ct.label))[sapply(core.en.regs, 
                                                                               "[[", "cell.type")])
qs::qsave(p.heatmap.fr.all, paste0(R.dir, "Heatmap_fraction_overlapped_enhs_all.qsave"))
save_image(p.heatmap.fr.all, paste0(image.dir, "Heatmap_fraction_overlapped_enhs_all.png"), 
           width = 2000, height = 1500)


p.heatmap.fr <- get_simple_heatmap(mat = core.dar.overlap[-(9:22), -(9:22)], splitted = NULL, name = "fr.overlap", 
                   bottom = 0.0, interium = 0.30, up = 1.0, labels = c("0.0", "1.0"),
                   at = c(0, 1), title_position = "leftcenter-rot", 
                   row.name.colors = setNames(col.list, 
                                              levels(act.obj$ct.label))[sapply(core.en.regs, 
                                                                               "[[", "cell.type")])
qs::qsave(p.heatmap.fr, paste0(R.dir, "Heatmap_fraction_overlapped_enhs.qsave"))
save_image(p.heatmap.fr, paste0(image.dir, "Heatmap_fraction_overlapped_enhs.png"), 
           width = 2000, height = 1500)



# Get B cell object
Idents(obj) <- obj$ct.label
b.obj <- subset(obj, ident = c("Normal B", "Tumor B prol", "Tumor B"))
dim(b.obj)
table(obj$ct.label)
sapply(cts.en.regs, "[[", "cell.type")
b.en.regs <- cts.en.regs[21:39]
sapply(b.en.regs, "[[", "cell.type")
b.de.genes <- FindAllMarkers(b.obj, assay = "RNA", 
                             features = Reduce("union", sapply(b.en.regs, "[[", "genes")))
b.sig.degs <- b.de.genes[b.de.genes$avg_log2FC > 0 & b.de.genes$p_val_adj < 0.05,]
dim(b.sig.degs)
b.da.enhs <- FindAllMarkers(b.obj, assay = "ATAC", min.pct = 0.05,
                             features = Reduce("union", sapply(b.en.regs, "[[", "enhancers")))
b.sig.dars <- b.da.enhs[b.da.enhs$avg_log2FC > 0 & b.da.enhs$p_val_adj < 0.05,]
dim(b.sig.dars)
qs::qsave(b.de.genes, paste0(R.dir, "DEGs_B_cells.qsave"))
qs::qsave(b.da.enhs, paste0(R.dir, "DARs_B_cells.qsave"))
b.en.regs <- lapply(b.en.regs, function(x) {
  x$de.genes <- intersect(x$genes, b.sig.degs[b.sig.degs$cluster == 
                                                x$cell.type,]$gene)
  x$da.enhs <- intersect(x$enhancers, b.sig.dars[b.sig.dars$cluster == 
                                                x$cell.type,]$gene)
  x
})
qs::qsave(b.en.regs, paste0(R.dir, "B_cell_eRegulons.qsave"))
sapply(b.en.regs, "[[", "cell.type")
sapply(b.en.regs, "[[", "de.genes") %>% sapply(., length)
sapply(b.en.regs, "[[", "da.enhs") %>% sapply(., length)
sapply(b.en.regs, "[[", "TF") %>% sapply(., length)
normal.en.regs <- b.en.regs[1:2]
prol.en.regs <- b.en.regs[3:10]
tumor.en.regs <- b.en.regs[11:19]



# Build eGRNs based on eRegulons on Tumor B prol
normal.en.grns <- eRegulons_to_eGRN(TFs = sapply(normal.en.regs, "[[", "TF"), 
                                  links = sapply(normal.en.regs, "[[", "links"), 
                                  col = "gene")
prol.en.grns <- eRegulons_to_eGRN(TFs = sapply(prol.en.regs, "[[", "TF"), 
                                  links = sapply(prol.en.regs, "[[", "links"), 
                                  col = "gene")
tumor.en.grns <- eRegulons_to_eGRN(TFs = sapply(tumor.en.regs, "[[", "TF"), 
                                  links = sapply(tumor.en.regs, "[[", "links"), 
                                  col = "gene")



# Calculate the weights of enh-gene links
b.links <- Reduce("union", pblapply(list(normal.en.grns, prol.en.grns, tumor.en.grns), 
                    function(x) {
                      paste0(x[, 2], "_", x[, 3]) %>% unique
                    }))
enh.gene.weights <- cal_enh_gene_across_ct(x = b.obj, y = c("Normal B", "Tumor B prol", "Tumor B"), 
                                           pairs = b.links, rna.assay = "RNA", atac.assay = "ATAC")
qs::qsave(enh.gene.weights, paste0(R.dir, "enh_gene_across_cells.qsave"))

b.tfEnhLinks <- Reduce("union", pblapply(list(normal.en.grns, prol.en.grns, tumor.en.grns), 
                                    function(x) {
                                      paste0(x[, 2], "_", x[, 1]) %>% unique
                                    }))
tf.enh.weights <- cal_enh_gene_across_ct(x = b.obj, y = c("Normal B", "Tumor B prol", "Tumor B"), 
                                         pairs = b.tfEnhLinks, rna.assay = "RNA", atac.assay = "ATAC")
qs::qsave(tf.enh.weights, paste0(R.dir, "TF_enh_across_cells.qsave"))


act.tf.enh.weights <- tf.enh.weights[apply(tf.enh.weights, 1, sum) > 0,]
rownames(act.tf.enh.weights) <- rownames(tf.enh.weights)[apply(tf.enh.weights, 1, sum) > 0]
dim(tf.enh.weights)
dim(act.tf.enh.weights)



# Colors
b.tfs <- c("TCF12", "SPI1", "STAT3", sapply(b.en.regs, "[[", "TF")) %>% unique
length(b.tfs)
color.ls <- setNames(
  c(
    "#D2001A", 
    "#367E18",
    "#5800FF",
    "#8758FF",
    "#FFB200",
    "#FA2FB5",
    "#7DCE13",
    "#66BFBF",
    "#5FD068",
    "#5EE6EB", 
    "#99A799",
    "#FFBC97"
  ), 
  c(b.tfs, "enhancer", "gene")
)
show_col(c(color.ls))


# Plotting for three Cytoscape images
rownames(enh.gene.weights)
rownames(act.tf.enh.weights)
cyto.tfs <- c("TCF12", "STAT3", "SPI1")
cyto.tf.enhs <- as.data.frame(act.tf.enh.weights)[strsplit(rownames(act.tf.enh.weights), 
                            split = "_") %>% 
                     sapply(., "[", 2) %in% cyto.tfs,]
rownames(cyto.tf.enhs) <- rownames(act.tf.enh.weights)[strsplit(rownames(act.tf.enh.weights), 
          split = "_") %>% 
    sapply(., "[", 2) %in% cyto.tfs]
cyto.enh.genes <- as.data.frame(enh.gene.weights)[strsplit(rownames(enh.gene.weights), 
                                            split = "_") %>% 
                                     sapply(., "[", 1) %in% 
                                     unique(strsplit(rownames(cyto.tf.enhs), split = "_") %>% 
                                     sapply(., "[[", 1)),]
rownames(cyto.enh.genes) <- rownames(enh.gene.weights)[strsplit(rownames(enh.gene.weights), 
                                                      split = "_") %>% 
                                               sapply(., "[", 1) %in% 
                                               unique(strsplit(rownames(cyto.tf.enhs), split = "_") %>% 
                                                        sapply(., "[[", 1))]
cyto.tf.genes <- union(strsplit(rownames(cyto.tf.enhs), split = "_") %>% sapply(., "[[", 2) %>% unique, 
                       strsplit(rownames(cyto.enh.genes), split = "_") %>% sapply(., "[[", 2) %>% unique)
cyto.enhs <- union(strsplit(rownames(cyto.tf.enhs), split = "_") %>% sapply(., "[[", 1) %>% unique, 
                   strsplit(rownames(cyto.enh.genes), split = "_") %>% sapply(., "[[", 1) %>% unique)
cyto.obj <- subset(b.obj, features = c(cyto.tf.genes, cyto.enhs))
length(cyto.tf.genes)
length(cyto.enhs)
dim(cyto.obj[['RNA']])
dim(cyto.obj[['ATAC']])



# Find variable features
DefaultAssay(cyto.obj) <- "RNA"
cyto.obj <- FindVariableFeatures(cyto.obj, selection.method = "vst", nfeatures = 300)
DefaultAssay(cyto.obj) <- "ATAC"
cyto.obj <- FindVariableFeatures(cyto.obj, selection.method = "vst", nfeatures = 300)
length(cyto.obj[['RNA']]@var.features)
length(cyto.obj[['ATAC']]@var.features)
var.enh.genes <- cyto.enh.genes[strsplit(rownames(cyto.enh.genes), split = "_") %>% 
                 sapply(., "[[", 1) %in% cyto.obj[['ATAC']]@var.features & 
                 strsplit(rownames(cyto.enh.genes), split = "_") %>% 
                 sapply(., "[[", 2) %in% cyto.obj[['RNA']]@var.features,]
var.tf.enhs <- cyto.tf.enhs[strsplit(rownames(cyto.tf.enhs), split = "_") %>% 
                                  sapply(., "[[", 1) %in% cyto.obj[['ATAC']]@var.features,]



# Co-binding
tfEnhDF <- Reduce("rbind", strsplit(rownames(var.tf.enhs), split = "_"))
colnames(tfEnhDF) <- c("enhancer", "TF")
spi1.stat3.tf.enh <- tfEnhDF[tfEnhDF[, 2] %in% c("SPI1", "STAT3"),]
spi1.tcf12.tf.enh <- tfEnhDF[tfEnhDF[, 2] %in% c("SPI1", "TCF12"),]
stat3.tcf12.tf.enh <- tfEnhDF[tfEnhDF[, 2] %in% c("STAT3", "TCF12"),]



# Co-regulation
enhGeneDF <- Reduce("rbind", strsplit(rownames(var.enh.genes), split = "_"))
colnames(enhGeneDF) <- c("enhancer", "gene")
spi1.stat3.enh.gene <- enhGeneDF[enhGeneDF[, 1] %in% spi1.stat3.tf.enh[, 1],]
spi1.tcf12.enh.gene <- enhGeneDF[enhGeneDF[, 1] %in% spi1.tcf12.tf.enh[, 1],]
stat3.tcf12.enh.gene <- enhGeneDF[enhGeneDF[, 1] %in% stat3.tcf12.tf.enh[, 1],]



# Functional genomics analysis for Tumor B Proliferation cells
cobind.gene.ls <- list(
  unique(spi1.stat3.enh.gene[, 2]), 
  unique(spi1.tcf12.enh.gene[, 2]),
  unique(stat3.tcf12.enh.gene[, 2])
)
func.ls <- run_GO_and_KEGG(genes.ll = cobind.gene.ls, dbs = "KEGG", 
                           org = "human")
sig.func <- pbmclapply(func.ls, function(x) {
  xx <- x[x$Adjusted.P.value < 0.05,]
  xx[order(xx$P.value),]
})
table(sig.func$KEGG_2019_Human$Id)
func.splitted <- split(sig.func$KEGG_2019_Human, sig.func$KEGG_2019_Human$Id) %>% sapply(., "[[", "Term")
Reduce("intersect", func.splitted)
setdiff(func.splitted[[1]], union(func.splitted[[2]], func.splitted[[3]]))
setdiff(func.splitted[[2]], union(func.splitted[[1]], func.splitted[[3]]))
setdiff(func.splitted[[3]], union(func.splitted[[1]], func.splitted[[2]]))
write.csv(sig.func$KEGG_2019_Human, paste0(table.dir, "Pathways_eGRN_SPI1_STAT3_TCF12.csv"), 
          quote = F, row.names = F)


spi1.stat3.pth <- c(
  "Leukocyte transendothelial migration",
  "Toxoplasmosis",
  "Human papillomavirus infection",
  "Thyroid hormone signaling pathway"
)


spi1.tcf12.pth <- c(
  "Acute myeloid leukemia",
  "Human cytomegalovirus infection",
  "Inositol phosphate metabolism",
  "Phospholipase D signaling pathway",
  "Choline metabolism in cancer"
)


stat3.tcf12.pth <- c(
  "Th1 and Th2 cell differentiation",
  "Inflammatory bowel disease (IBD)"
)



pth.table <- sig.func$KEGG_2019_Human
colnames(pth.table)
pth.table$Id
spi1.stat3.table <- pth.table[pth.table$Id == "1",]
spi1.tcf12.table <- pth.table[pth.table$Id == "2",]
stat3.tcf12.table <- pth.table[pth.table$Id == "3",]
dim(spi1.stat3.table)
dim(spi1.tcf12.table)
dim(stat3.tcf12.table)

spi1.stat3.table[1:5, 1:5]
spi1.tcf12.table[1:5, 1:5]
stat3.tcf12.table[1:5, 1:5]

spi1.stat3.key.table <- spi1.stat3.table[spi1.stat3.table$Term %in% spi1.stat3.pth,]
spi1.tcf12.key.table <- spi1.tcf12.table[spi1.tcf12.table$Term %in% spi1.tcf12.pth,]
stat3.tcf12.key.table <- stat3.tcf12.table[stat3.tcf12.table$Term %in% stat3.tcf12.pth,]


# Select example DEGs in Tumor B proliferation cells and
# enriched in important pathways
de.genes
range(de.genes$p_val_adj)
table(de.genes$cluster)
normal.degs <- de.genes[de.genes$cluster == "Normal B",]
prol.degs <- de.genes[de.genes$cluster == "Tumor B prol",]
tumor.degs <- de.genes[de.genes$cluster == "Tumor B",]
dim(de.genes)
dim(normal.degs)
dim(prol.degs)
dim(tumor.degs)



# Normal B
normal.gD <- df_to_Cyto(T2R = var.tf.enhs[, 1, drop = F], 
                        R2G = var.enh.genes[, 1, drop = F])
vcount(normal.gD)
ecount(normal.gD)
plot(normal.gD)



# Tumor B prol
prol.gD <- df_to_Cyto(T2R = var.tf.enhs[, 2, drop = F], 
                        R2G = var.enh.genes[, 2, drop = F])
vcount(prol.gD)
ecount(prol.gD)
plot(prol.gD)



# Tumor B
tumor.gD <- df_to_Cyto(T2R = var.tf.enhs[, 3, drop = F], 
                      R2G = var.enh.genes[, 3, drop = F])
vcount(tumor.gD)
ecount(tumor.gD)
plot(tumor.gD)
qs::qsave(list(colors = color.ls, normal = normal.gD, 
               prol = prol.gD, tumor = tumor.gD), paste0(R.dir, "Vis_igraphs.qsave"))



# Exploit differential regulation in Tumor B Proliferation cells
spi1.stat3.pth.genes <- strsplit(spi1.stat3.key.table$Genes, split = ";")
spi1.tcf12.pth.genes <- strsplit(spi1.tcf12.key.table$Genes, split = ";")
stat3.tcf12.pth.genes <- strsplit(stat3.tcf12.key.table$Genes, split = ";")

spi1.stat3.rg <- Reduce("rbind", pblapply(spi1.stat3.pth.genes, function(x) {
  c(
    strength(normal.gD, x, mode = "all") %>% mean, 
    strength(prol.gD, x, mode = "all") %>% mean, 
    strength(tumor.gD, x, mode = "all") %>% mean
  )
}))
spi1.stat3.rg

spi1.tcf12.rg <- Reduce("rbind", pblapply(spi1.tcf12.pth.genes, function(x) {
  c(
    strength(normal.gD, x, mode = "all") %>% mean, 
    strength(prol.gD, x, mode = "all") %>% mean, 
    strength(tumor.gD, x, mode = "all") %>% mean
  )
}))
spi1.tcf12.rg

stat3.tcf12.rg <- Reduce("rbind", pblapply(stat3.tcf12.pth.genes, function(x) {
  c(
    strength(normal.gD, x, mode = "all") %>% mean, 
    strength(prol.gD, x, mode = "all") %>% mean, 
    strength(tumor.gD, x, mode = "all") %>% mean
  )
}))
stat3.tcf12.rg


# Make the table using the top-30 rows
require(DescTools)
dim(var.tf.enhs)
dim(var.enh.genes)
head(var.tf.enhs)
head(var.enh.genes)
write.csv(var.enh.genes, quote = F, paste0(table.dir, "enh_gene_scores.csv"))
var.tf.enhs.named <- var.tf.enhs
rownames(var.tf.enhs.named) <- paste0(
  Reduce("rbind", strsplit(rownames(var.tf.enhs), split = "_"))[, 2], "_", 
  Reduce("rbind", strsplit(rownames(var.tf.enhs), split = "_"))[, 1])
write.csv(var.tf.enhs.named, 
  quote = F, paste0(table.dir, "TF_enh_scores.csv"))


save.image(paste0(scratch.dir, "5_DSLL_global_analyses.RData"))
# save.image("/fs/ess/scratch/PCON0022/liyang/stream/Case_1_DSLL_cell_type/Rdata/5_DSLL_global_analyses.RData")
