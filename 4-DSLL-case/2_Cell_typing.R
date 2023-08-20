########################################################################
#                                                                      #
# 1. UMAP: cell types in which eRegulons are active                    #
# 2. Piechart: eRegulon numbers and percentages in various cell types  #
# 3. Barplots: RSS scores of eRegulons, genes, enhancers, and linkages #
# 4. Dotted heatmap: TF expression - RSS scores                        #
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
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
setwd(R.dir)


# Set colors
full.colors <- c(
  "#FF0000", 
  "#FEB139",
  "#FAD9A1",
  "#66BFBF",
  "#BD4291",
  "#187498",
  "#143F6B",
  "#357C3C",
  "#FF0075",
  "#30475E",
  "#FFB400", 
  "#3AB4F2",
  "#764AF1",
  "#F66B0E",
  "#066163",
  "#B20600", 
  "#8D8DAA", 
  "#B85252", 
  "#6C4A4A",
  "#16A596", 
  "#FA7F72", 
  "#EBEBEB", 
  "#F7DAD9", 
  "#19D3DA", 
  "#32E0C4", 
  "#E5DF88", 
  "#FFCB8E", 
  "#D49A89", 
  "#CEDEBD",
  "#99B19C"
)




# Load data
en.regs <- qs::qread(paste0(R.dir, "Extended_eRegs_coh_3.qsave"))
names(en.regs[[1]])
sapply(en.regs, "[[", "links") %>% sapply(., length) %>% range # 15 1246
act.cells <- Reduce("union", sapply(en.regs, "[[", "cells"))
length(act.cells) # 3846
obj <- qs::qread(paste0(R.dir, "lymph_obj_annotated.qsave"))
dim(obj) # 8770
intersect(act.cells, colnames(obj)) %>% length # 3846


# Get the subsetted object on active cells
act.obj <- subset(x = obj, cells = act.cells)
dim(act.obj) # 3846
qs::qsave(act.obj, paste0(R.dir, "Obj_active.qsave"))
table(act.obj$cell.type)
table(act.obj$ct.label)




# Identify cell-type enrichment for eRegulons against cell types from all cells
act.cell.list <- sapply(en.regs, "[[", "cells")
names(act.cell.list) <- seq_along(act.cell.list)
celltypes <- readRDS(paste0(R.dir, "DSLL_DeepMAPS_celltypes.RDS"))
celltype.list <- split(celltypes, f = celltypes) %>% sapply(., names)
names(celltype.list) <- levels(obj$ct.label)
hyper.df <- run_hyper_test(x = act.cell.list, y = celltype.list, n = length(celltypes))
sig.pairs <- hyper.df[hyper.df$Pval * length(celltype.list) < 0.05,]
dim(sig.pairs) # 39
unique(sig.pairs$x) %>% length # 23
names(celltype.list)[unique(sig.pairs$y)]




# Identify cell-type enrichment for eRegulons against cell types from active cells
act.celltypes <- celltypes[act.cells]
act.celltype.list <- split(act.celltypes, f = act.celltypes) %>% sapply(., names)
hyper.df <- run_hyper_test(x = act.cell.list, y = act.celltype.list, n = length(act.celltypes))
sig.pairs <- hyper.df[hyper.df$Pval * length(celltype.list) < 0.05,]
sig.pairs <- cbind(sig.pairs, padj = sig.pairs$Pval * length(celltype.list))
sig.pairs$y <- names(celltype.list)[sig.pairs$y]
dim(sig.pairs)
unique(sig.pairs$x) %>% length 
b.celltypes <- c("Normal B", "Tumor B prol", "Tumor B")
sig.pairs[sig.pairs$y %in% b.celltypes, c(1, 2)] %>% split(., .$y) %>% sapply(., nrow)
sig.pairs[sig.pairs$y == "Normal B", "x"] 
sig.pairs[sig.pairs$y == "Tumor B prol", "x"] 
sig.pairs[sig.pairs$y == "Tumor B", "x"] 
sapply(en.regs, "[[", "TF")[sig.pairs[sig.pairs$y == "Normal B", "x"]] 
sapply(en.regs, "[[", "TF")[sig.pairs[sig.pairs$y == "Tumor B prol", "x"]]
sapply(en.regs, "[[", "TF")[sig.pairs[sig.pairs$y == "Tumor B", "x"]]


# Get cell-type-specific eGRNs
ct.en.reg.list <- split(sig.pairs[, 1:2], f = sig.pairs[, 2]) %>% sapply(., "[[", "x")
length(ct.en.reg.list) 
names(ct.en.reg.list)
cts.en.regs <- do.call("c", pbmclapply(seq_along(ct.en.reg.list), function(i) {
  x <- ct.en.reg.list[[i]]
  TF.id <- data.frame(TF = sapply(en.regs, "[[", "TF")[x], 
                      id = seq_along(sapply(en.regs, "[[", "TF")[x])) %>% split(., .$TF) %>% 
    sapply(., "[[", "id")
  lapply(seq_along(TF.id), function(j) {
    y <- TF.id[[j]]
    TF <- names(TF.id)[j]
    genes <- Reduce("union", sapply(en.regs[x][y], "[[", "genes"))
    enhancers <- Reduce("union", sapply(en.regs[x][y], "[[", "peaks"))
    links <- merge_GRanges(sapply(en.regs[x][y], "[[", "links"), 
                           col = "gene")
    cell.type <- names(ct.en.reg.list)[i]
    list(TF = TF, genes = genes, enhancers = enhancers, links = links, 
         cell.type = cell.type)
  })
}, mc.cores = detectCores()))
sapply(cts.en.regs, "[[", "TF")
sapply(cts.en.regs, "[[", "cell.type")
sapply(en.regs, "[[", "TF")
qs::qsave(cts.en.regs, paste0(R.dir, "Cell_type_specific_eRegulons.qsave"))




# Compute TF expression across all cells
rna.m <- GetAssayData(obj, slot = "data", assay = "RNA")
dim(rna.m)
cts.en.regs <- pbmclapply(cts.en.regs, function(x) {
  if (!x$TF %in% rownames(rna.m)) {
    x$TF.exp <- 0
  } else {
    x$TF.exp <- rna.m[x$TF, rownames(obj@meta.data[obj@meta.data$ct.label == x$cell.type,])] %>% mean
  }
  return(x)
}, mc.cores = detectCores())
sapply(cts.en.regs, "[[", "TF.exp") %>% range 
sapply(cts.en.regs, "[[", "cell.type")




# Compute TF expression across active cells
act.rna.m <- GetAssayData(act.obj, slot = "data", assay = "RNA")
dim(act.rna.m) 
cts.en.regs <- pbmclapply(cts.en.regs, function(x) {
  if (!x$TF %in% rownames(act.rna.m)) {
    x$TF.act.exp <- 0
  } else {
    x$TF.act.exp <- act.rna.m[x$TF, rownames(act.obj@meta.data[act.obj@meta.data$ct.label == 
                                                                 x$cell.type,])] %>% mean
  }
  return(x)
}, mc.cores = detectCores())
sapply(cts.en.regs, "[[", "TF.act.exp") %>% range 
sapply(cts.en.regs, "[[", "cell.type")


# Get TF expression matrix (TF x cell type)
TF.exp <- Reduce("rbind", pbmclapply(cts.en.regs, function(x) {
  sapply(unique(sapply(cts.en.regs, "[[", "cell.type")), function(y) {
    rna.m[x$TF, rownames(obj@meta.data[obj@meta.data$ct.label == y,])] %>% mean
  })
}, mc.cores = detectCores()))
rownames(TF.exp) <- seq_along(cts.en.regs)
TF.exp[1:3, 1:3]
heatmap(TF.exp, scale = "row", Colv = NA, Rowv = NA)
heatmap(TF.exp[c(21, 22, 30:46), c(6, 8, 9)], scale = "row", Colv = NA, Rowv = NA)


# DEG and DAR analyses between cell types on all cells
Idents(act.obj) <- act.obj$ct.label
de.genes <- FindAllMarkers(act.obj, assay = "RNA", 
                           features = union(Reduce("union", sapply(cts.en.regs, "[[", "genes")), 
                                            Reduce("union", sapply(cts.en.regs, "[[", "TF"))))
dim(de.genes)    
da.enhs <- FindAllMarkers(act.obj, assay = "ATAC", min.pct = 0.01,
                           features = Reduce("union", sapply(cts.en.regs, "[[", "enhancers")))
dim(da.enhs)          
qs::qsave(de.genes, paste0(R.dir, "DEGs.qsave"))
qs::qsave(da.enhs, paste0(R.dir, "DARs.qsave"))




# Add DEGs and DARs to cell-type-specific eRegulons
de.genes <- de.genes[de.genes$avg_log2FC > 0.25 &
                       de.genes$p_val_adj < 0.05,]
dim(de.genes) 
da.enhs <- da.enhs[da.enhs$avg_log2FC > 0.00 &
                     da.enhs$p_val_adj < 0.05,]
dim(da.enhs)     
de.genes <- de.genes[with(de.genes, order(p_val, -avg_log2FC)),]
da.enhs <- da.enhs[with(da.enhs, order(p_val, -avg_log2FC)),]
head(de.genes)
head(da.enhs)
cts.en.regs <- pbmclapply(cts.en.regs, function(x) {
  gene.df <- de.genes[de.genes$cluster == x$cell.type & de.genes$gene %in% x$genes,]
  enh.df <- da.enhs[da.enhs$cluster == x$cell.type & da.enhs$gene %in% x$enhancers,]
  x$de.genes <- head(gene.df$gene, n = Inf)
  x$da.enhs <- head(enh.df$gene, n = Inf)
  if (x$TF %in% de.genes[de.genes$cluster == x$cell.type & de.genes$p_val_adj < 0.05 & 
                         de.genes$avg_log2FC > 0.25,]$gene) {
    x$de.TF <- T
  } else {
    x$de.TF <- F
  }
  x
}, mc.cores = detectCores())
lapply(cts.en.regs, "[[", "de.genes") %>% sapply(., length)
lapply(cts.en.regs, "[[", "da.enhs") %>% sapply(., length)
sapply(cts.en.regs, "[[", "cell.type")
sapply(cts.en.regs, "[[", "TF")[sapply(cts.en.regs, "[[", "de.TF")] # "TCF7"  "SPI1"  "TCF12"




# Build heatmap for DEGs in eRegulons
act.obj <- ScaleData(act.obj)
exp.m <- act.obj@assays$RNA@data
top10.exp <- Reduce("rbind", pbmclapply(cts.en.regs, function(x) {
  sapply(levels(act.obj$ct.label), function(y) {
    mean(exp.m[x$de.genes, names(act.obj$ct.label)[act.obj$ct.label == y]])
  })
}))
rownames(top10.exp) <- seq_along(cts.en.regs)
heatmap(top10.exp, scale = "row", Rowv = NA, Colv = NA)
sapply(cts.en.regs, "[[", "cell.type")




# Plot 1: generate the UMAP plots for all cells and eRegulon-active cells


get_UMAP <- function(object = NULL, reduction.method = "wnn.umap", 
                     pt_size = 0.5, txt = "Phenotype",
                     combined.cell.types = NULL, 
                     legend.position = "right", 
                     path = NULL, color.ll = NULL,
                     margin.quad = c(10.0, -5.0, 0.0, 0.0), 
                     ...) {
  
  # Check whether embedding exists
  if (!"Seurat" %in% class(object)) {
    stop ("This object is not a Seurat object!\n", 
          "Please check the method to generate this object!\n")
  }
  meta.data <- object[[]]
  message ("This Seurat object has the following meta data: \n", 
           "====================================\n", 
           paste(colnames(meta.data), collapse = ", "), 
           "\n====================================\n")
  if (!reduction.method %in% names(object@reductions)) {
    warning ("This Seurat object does not contain the result of reduction method ", 
             reduction.method, ".\n", 
             "We will generate the reduction result now ...\n")
    object <- Seurat_reduce_dim(object = object, reduction.method = reduction.method, 
                                n.dim = 50) # dimension reduction via Seurat
  }
  
  
  # Save the number of cells in each cell type
  my.plot.all.source <- cbind.data.frame(Embeddings(object, reduction = reduction.method),
                                         Cell_type = Idents(object))
  tmp.celltype <- levels(unique(my.plot.all.source$Cell_type)) # all phenotypes
  
  
  # Generate the image
  require(ggplot2)
  p.cluster <- ggplot(my.plot.all.source, 
                      aes(x = my.plot.all.source[, 1], y = my.plot.all.source[, 2])) + 
    xlab(colnames(my.plot.all.source)[1]) + ylab(colnames(my.plot.all.source)[2]) # an empty image
  message ("Generating the coordinates for plotting ...\n")
  p.cluster <- p.cluster + geom_point(stroke = pt_size, 
                                      size = pt_size, aes(col = my.plot.all.source[, "Cell_type"]))
  message ("Add points representing cells to the coordinates ...\n")
  # add points and legends
  
  p.cluster <- p.cluster + guides(colour = guide_legend(override.aes = list(size = 5))) # no change
  
  
  # Change the colors of points
  message ("Changing the colors of points ...\n")
  color.ll <- full.colors
  # if (!is.null(color.ll)) {
  #   color.ll <- full.colors
  # } else {
  #   require(Polychrome)
  #   color.ll <- as.character(palette36.colors(36)[-2][1 : length(tmp.celltype)])
  # } else if (length(tmp.celltype) < 5) {
  #   p.cluster <- p.cluster + 
  #     scale_colour_manual(name = paste(txt, ":(Cells)", sep = ""), 
  #                         values = brewer.pal(4, "Spectral")[c(2, 1, 3, 4)],
  #                         breaks = tmp.celltype,
  #                         labels = paste0(tmp.celltype,":(", 
  #                                         as.character(summary(my.plot.all.source$Cell_type)), ")"))
  # } else {
  #   if (length(tmp.celltype) == 2) {
  #     p.cluster <- p.cluster + 
  #       scale_colour_manual(name = paste(txt, ": (Cells)", sep = ""), 
  #                           values = color.ll[c(1, length(tmp.celltype))], 
  #                           breaks = tmp.celltype, 
  #                           labels = paste0(tmp.celltype, ": (", 
  #                                           as.character(summary(my.plot.all.source$Cell_type)), 
  #                                           ")"))
  #   } else if (length(tmp.celltype) == 3) {
  #     p.cluster <- p.cluster + 
  #       scale_colour_manual(name = paste(txt, ": (Cells)", sep = ""), 
  #                           values = color.ll, 
  #                           breaks = tmp.celltype, 
  #                           labels = paste0(tmp.celltype, ": (", 
  #                                           as.character(summary(my.plot.all.source$Cell_type)), 
  #                                           ")"))
  #   }
  # }
  
  
  # Post processing
  require(formattable)
  p.cluster <- p.cluster + 
    scale_colour_manual(name = paste(txt, ": (Cells)", sep = ""), 
                        values = color.ll,
                        breaks = tmp.celltype, 
                        labels = paste0(tmp.celltype, ": (", 
                                        as.character(comma(summary(my.plot.all.source$Cell_type), 
                                                           format = "d")), ")"))
  
  message ("The following colors are used to generate the UMAP: \n", 
           "=================================================\n", 
           paste(color.ll, collapse = "\n"), 
           "\n=================================================\n\n", 
           paste0(tmp.celltype, collapse = "\n"), 
           "\n\n")
  p.cluster <- p.cluster + theme_classic()
  p.cluster <- p.cluster + coord_fixed(ratio = 1)
  p.cluster <- p.cluster + theme(legend.position = legend.position)
  p.cluster <- p.cluster + theme(plot.margin = unit(margin.quad, "pt"))
  save_image(p = p.cluster, path = path)
}


Idents(obj) <- obj$ct.label
p.umap <- get_UMAP(object = obj, reduction.method = "umap.rna", txt = "Cell type")
qs::qsave(p.umap, paste0(R.dir, "UMAP_DeepMAPS_celltypes.qsave"))
Idents(act.obj) <- act.obj$ct.label
p.act.umap <- get_UMAP(object = act.obj, reduction.method = "umap.rna", txt = "Cell type")
qs::qsave(p.act.umap, paste0(R.dir, "UMAP_DeepMAPS_active_celltypes.qsave"))




# Plot 2: generate the barplot of eRegulons enriched in cell types

get_bar_plot <- function(df = NULL, fill = "#757575", 
                         variable = variable, label = value, fill.name = NULL,
                         value = value, group = group, position = "dodge", x.lab = "x", 
                         y.lab = "y", 
                         angle = 0, size = 10, path = NULL, width = 500, height = 500) {
  
  require(ggplot2)
  p <- ggplot(df, aes(x = variable, y = value, fill = as.factor(group))) + 
    geom_bar(stat = "identity", position = position) + xlab(label = x.lab) + 
    ylab(label = y.lab) + geom_text(aes(label = sprintf("%0.2f", value) ), vjust = -0.50) + 
    theme(axis.text.x = element_text(angle = angle, colour = "black", family = "Arial", size = size), 
          axis.text.y = element_text(colour = "black", family = "Arial", size = size), 
          axis.title.x = element_text(colour = "black", family = "Arial", size = 1.2 * size), 
          axis.title.y = element_text(colour = "black", family = "Arial", size = 1.2 * size), 
          panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.background = element_blank(), legend.position = "none",
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_fill_manual(values = fill, name = fill.name)
  save_image(p = p, path = path, width = width, height = height)
}


ct.eReg.df <- sapply(split(sig.pairs[, 1:2], f = sig.pairs[, 2]), nrow) %>% as.data.frame
ct.eReg.df <- cbind("Cell type" = rownames(ct.eReg.df), ct.eReg.df)
rownames(ct.eReg.df) <- NULL
colnames(ct.eReg.df) <- c("Cell type", "eRegulon")
ct.eReg.df <- ct.eReg.df %>%
  arrange(sapply(`Cell type`, function(y) which(y == levels(act.obj$ct.label)[1:9])))
ct.eReg.df$`Cell type` <- factor(ct.eReg.df$`Cell type`, levels = levels(act.obj$ct.label)[1:9])
p.bar.ct.eRegs <- get_bar_plot(df = ct.eReg.df, 
                               fill = full.colors[1:length(levels(act.obj$ct.label)[1:9])])
qs::qsave(p.bar.ct.eRegs, paste0(R.dir, "Barplot_celltype_eRegulons.qsave"))




save.image("2_DSLL_celltype_eRegulons.RData")
