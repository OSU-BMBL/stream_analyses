#####################################################
#                                                   #
#                   EX analysis                     #
#                                                   #
#####################################################



# Libraries
require(easypackages)
libs <- c(
  "Seurat", 
  "Signac", 
  "dplyr", 
  "parallel", 
  "pbmcapply", 
  "pbapply",
  "scales", 
  "data.table", 
  "monocle3"
)
libraries(libs)



# Parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_3_AD/Rdata/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_3_AD/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_3_AD/Images/"
scratch.dir <- "/fs/ess/scratch/PCON0022/liyang/stream/Case_2_AD/"
setwd(R.dir)



###############################################################
#                                                             #
#  4. monocle3_TI_atac : trajectory inference for ATAC data   #
#                                                             #
###############################################################


monocle3_TI_atac <- function(mat = NULL, root_cells = NULL, 
                             pseudotime.annot = NULL, nBins = 20, 
                             binarize = F) {
  
  # Load data
  message ("Loading the accessibility data: ")
  dim(mat)
  
  
  # Convert matrix into data.frame
  message ("Converting matrix into triplets ...")
  require(monocle3)
  require(cicero)
  if (grepl("-", rownames(mat)[1])) {
    message ("Substituting '-' with '_' for row names ...")
    rownames(mat) <- gsub("-", "_", rownames(mat))
  }
  # if (grepl("-", colnames(mat)[1])) {
  #   message ("Substituting '-' with '_' for column names ...")
  #   colnames(mat) <- gsub("-", "_", colnames(mat))
  # }
  cicero_data <- as.data.frame(summary(mat))
  colnames(cicero_data) <- paste0("V", 1:3)
  cicero_data$V1 <- rownames(mat)[cicero_data$V1]
  cicero_data$V2 <- colnames(mat)[cicero_data$V2]
  input_cds <- make_atac_cds(cicero_data)
  
  
  # Order cells
  set.seed(2017)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                                preprocess_method = "LSI")
  input_cds <- cluster_cells(input_cds)
  input_cds <- learn_graph(input_cds)
  if (!is.null(pseudotime.annot)) {
    if (grepl("-", names(pseudotime.annot)[1])) {
      message ("Adding pseudotime information ...")
      names(pseudotime.annot) <- gsub("-", "_", names(pseudotime.annot))
    }
    pData(input_cds)$Pseudotime <- pseudotime.annot
  }else if (is.null(root_cells)) {
    message ("Inferring trajectory using annotated root cells ...")
    input_cds <- order_cells(input_cds, 
                             root_cells = root_cells)
  } else {
    message ("Inferring trajectory without annotated root cells ...")
    input_cds <- order_cells(input_cds)
  }
  
  
  # First, assign a column in the pData table to umap pseudotime
  message ("Assigning cells to bins by cutting the pseudotime trajectory into ", 
           nBins, " bins ...")
  input_cds_lin <- input_cds
  # pData(input_cds_lin)$Pseudotime <- pseudotime(input_cds_lin)
  pData(input_cds_lin)$cell_subtype <- cut(pData(input_cds_lin)$Pseudotime, nBins)
  # pData(input_cds_lin)$cell_subtype <- cut(pseudotime(input_cds_lin), bin.size)
  binned_input_lin <- aggregate_by_cell_bin(input_cds_lin, "cell_subtype")
  
  
  # For speed, run fit_models on 1000 randomly chosen genes
  message ("Identifying peaks whose accessibilities change as a function of pseudotime ...")
  set.seed(1000)
  acc_fits <- fit_models(binned_input_lin, 
                         model_formula_str = "~Pseudotime + num_genes_expressed" )
  fit_coefs <- coefficient_table(acc_fits)
  message ("There are in total ", length(which(fit_coefs$q_value < 0.05)), 
           " peaks significant under q-value: 0.05.\n", 
           "There are in total ", length(which(fit_coefs$p_value < 0.05)), 
           " peaks significant under p-value: 0.05.\n")
  
  
  # Change names
  if (grepl("_", rownames(input_cds_lin)[1])) {
    message ("Substituting '_' with '-' ...")
    rownames(input_cds_lin) <- gsub("_", "-", rownames(input_cds_lin))
    fit_coefs$site_name <- gsub("_", "-", fit_coefs$site_name)
  }
  
  
  list(cds = input_cds_lin, da.regions = fit_coefs)
}


save_image <- function(p, width = 1500, height = 1000, path = NULL) {
  
  
  # p : an object in the formats of gg, ggplot, or ggarrange
  # width : the width of the image, 1500 as default
  # height : the height of the image, 1000 by default
  # path : the path to save the image
  # text.font : the font in figure, e.g., text.font <- "Arial"
  
  
  # Check the parameters
  arrange.formats <- c("egg", "gtable", "gTree", "grob", "gDesc")
  if (length(intersect(class(p), arrange.formats)) > 0) {
    message ("Converting the format from ", class(p)[1], " into ggplot.\n")
    library(ggplot2)
    p <- as.ggplot(p)
  }
  if (length(intersect(class(p), c("gg", "ggplot", "ggarrange", 
                                   "Heatmap", "ComplexHeatmap", 
                                   "ggsurvplot"))) < 1) {
    stop ("The input object is not an image in the formats of gg, ggplot, ggarrange, Heatmap, and ComplexHeatmap!\n", 
          "Please check the input object!.\n")
  }
  
  
  # Print images
  if (is.null(path)) {
    message ("Finished saving the image in the memory.\n")
    if (!"ggplot"%in% class(p) & 
        class(p)[1] != "ggsurvplot") {
      library(ggplotify)
      p <- as.ggplot(p)
    }
    return(p)
  } else if (grepl(".eps$", path, ignore.case = T)) {
    postscript(path, width = width, 
               height = height)
    print(p)
    dev.off()
  } else if (grepl(".tiff$", path, ignore.case = T)) {
    tiff(filename = path, width = width, height = height, res = 300)
    print(p)
    dev.off()
  } else if (grepl(".png$", path, ignore.case = T)) {
    png(filename = path, width = width, height = height, res = 300)
    print(p)
    dev.off()
  }
  message ("The image has been saved to ", path, ".\n")
}


get_simple_heatmap <- function(mat = NULL, splitted = NULL, name = NULL, rowfont = 20,
                               path = NULL, low = "#21209C", border.col = "black", 
                               middle = "#FF5F00", labels = c("0.0", "1.0"),
                               at = c(0, 1), title_position = "leftcenter-rot",
                               high = "yellow", bottom = 0.0, interium = 0.5, 
                               row.name.colors = NULL, legend.fontsize = 20,
                               up = 1.0) {
  
  # mat : the input data frame
  # splitted : the vector representing different groups of rows
  # name : the name of legends
  # path : the path to save the image
  # low : the color representing low level
  # high : the color denoting high level
  
  message ("The matrix contains ", nrow(mat), " rows and ", 
           ncol(mat), " columns.\n")
  require(ComplexHeatmap)
  require(circlize)
  col_fun = colorRamp2(c(bottom, interium, up), c(low, middle, high))
  if (!is.null(splitted)) {
    p <- Heatmap(mat, row_split = splitted, show_row_dend = F, column_labels = F,
                 row_title = NULL, name = name, cluster_rows = F,
                 show_column_dend = F, cluster_columns = F, col = col_fun, 
                 rect_gp = gpar(col = NA, lwd = 0.5), show_column_names = F, 
                 heatmap_legend_param = list(
                   title_position = title_position,
                   at = at, 
                   border = border.col, 
                   labels = labels
                 ),
                 row_names_gp = gpar(col = if (is.null(row.name.colors)) {
                   rep("black", nrow(mat))
                 }else {
                   row.name.colors
                 }),
                 left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = full.colors),
                                                                  labels = unique(splitted),
                                                                  labels_gp = gpar(col = "white",
                                                                                   fontsize = 10))))
  } else {
    p <- Heatmap(mat, show_row_dend = F, cluster_rows = F, 
                 row_names_side = "left", show_column_names = F, 
                 show_column_dend = F, cluster_columns = F, 
                 heatmap_legend_param = list(
                   title_position = title_position,
                   at = at, 
                   border = border.col, 
                   labels = labels, 
                   title_gp = gpar(fontsize = legend.fontsize * 1.0),
                   labels_gp = gpar(fontsize = legend.fontsize)
                 ), 
                 row_names_gp = gpar(col = if (is.null(row.name.colors)) {
                   rep("black", nrow(mat))
                 } else {
                   row.name.colors
                 }, 
                 fontsize = rowfont),
                 row_title = NULL, name = name, 
                 col = col_fun)
  }
  
  save_image(p = p, path = path)
}


###############################################################
#                                                             #
# 29. jaccard_pseudotime : generate dotplot of                #
#     Jaccard-pseudotime plot                                 #
#                                                             #
###############################################################


# Input :
# 1. x : the first variable, e.g., expression
# 2. y : the second variable, e.g., accessibility
# 3. z : the third variable, e.g., pseudotime
# 4. nBins : the number of bins
# 5. div.method : the method to calculate Jaccard index
# 6. path : the path tyo save the image


# Generate the fitting formula
eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}


# Main function
jaccard_pseudotime <- function(x = NULL, y = NULL, fontsize = 20,
                               z = NULL, nBins = 50, div.method = "min", ifPlot = T,
                               path = NULL, eq.y = 0.3, eq.x = 1.5, rsquare.x = 1.5,
                               rsquare.y = 0.35) {
  
  require(dplyr)
  message ("Cutting pseudotime into ", nBins, " bins ...")
  df <- data.frame(cells = names(z), bin = cut(z, nBins)) %>% 
    split(., .$bin) %>% sapply(., "[[", "cells")
  
  
  # Calculating Jaccard indices
  require(pbmcapply)
  res <- pbmclapply(df, function(i) {
    xx <- x[i]
    yy <- y[i]
    xxx <- names(which(xx > 0))
    yyy <- names(which(yy > 0))
    if (div.method == "union") {
      div <- length(union(xxx, yyy))
    } else if (div.method == "min") {
      div <- min(length(xxx), length(yyy))
    } else if (div.method == "max") {
      div <- max(length(xxx), length(yyy))
    }
    length(intersect(xxx, yyy)) / div
  }, mc.cores = detectCores()) %>% unlist
  jdx <- res
  jdx[is.na(jdx)] <- 0
  pt <- strsplit(names(jdx), split = ",") %>% sapply(., "[[", 2) %>% gsub(".$", "", .) %>% 
    as.numeric
  
  
  # Plotting
  require(ggplot2)
  require(ggpubr)
  corr <- cor.test(x = pt, y = jdx, method = 'pearson')
  # if (pos.neg == 'pos' & corr$estimate < 0 | pos.neg == 'neg' & corr$estimate > 0) {
  #   return(0)
  # }
  if (!ifPlot) {
    message ('We will not generate plotting ...')
    return(list(jdx, corr$p.value, corr$estimate))
  }
  data <- data.frame(Pseudotime = pt, Jaccard = jdx)
  p <- ggplot(data, aes(x = Pseudotime, y = Jaccard)) +
    geom_point() +
    geom_smooth(method = lm , color = "red", fill = "#69b3a2", se = TRUE) + 
    ylab("Fraction of cells") +  stat_regline_equation(label.y = eq.y, label.x = eq.x, 
                                                       aes(label = ..eq.label..), size = round(fontsize * .3)) +
    stat_regline_equation(label.y = rsquare.y, label.x = rsquare.x, 
                          aes(label = ..rr.label..), size = round(fontsize * .3)) + 
    theme_bw() + theme(axis.text = element_text(size = fontsize, colour = 'black'), 
                       axis.title = element_text(size = fontsize, colour = 'black'))
  # geom_text(x = axis.x, y = axis.y, label = eq(data$Pseudotime, data$Jaccard), parse = T)
  
  
  save_image(p = p, path = path)
}



# Parameters
ct <- "IN"



# Load data
obj <- qs::qread(paste0(R.dir, "Object_Shane_annot.qsave"))
Idents(obj) <- obj$month
table(obj$month)
obj.ls <- qs::qread(paste0(R.dir, "obj_list.qsave"))
sapply(obj.ls, ncol)
obj.2.5 <- obj.ls$`2.5`
obj.5.7 <- obj.ls$`5.7`
obj.13 <- obj.ls$`13`
en.reg.2.5 <- qs::qread(paste0(R.dir, "Extended_eGRNs_2.5.qsave"))
en.reg.5.7 <- qs::qread(paste0(R.dir, "Extended_eGRNs_5.7.qsave"))
en.reg.13 <- qs::qread(paste0(R.dir, "Extended_eGRNs_13.qsave"))
sig.cts.2.5 <- qs::qread(paste0(R.dir, "sig_cts_2.5.qsave"))
sig.cts.5.7 <- qs::qread(paste0(R.dir, "sig_cts_5.7.qsave"))
sig.cts.13 <- qs::qread(paste0(R.dir, "sig_cts_13.qsave"))



# Get cell-type-specific eRegulons in EX
sig.cts.2.5
sig.cts.5.7
sig.cts.13
cts.en.regs <- c(
  en.reg.2.5[unique(sig.cts.2.5[sig.cts.2.5$y == ct, "x"])],
  en.reg.5.7[unique(sig.cts.5.7[sig.cts.5.7$y == ct, "x"])],
  en.reg.13[unique(sig.cts.13[sig.cts.13$y == ct, "x"])]
)
length(cts.en.regs)


# Calculate fraction of nonzero expression and accessibility submartrices
links <- pbmclapply(cts.en.regs, function(x) {
  paste0(GRangesToString(x$links), "_", x$links$gene)
}, mc.cores = detectCores())
length(links)
sapply(links, length)
u.links <- Reduce("union", links)
length(u.links)


# Calculate link matrix
link.mat <- Reduce("rbind", pbmclapply(u.links, function(x) {
  ar <- strsplit(x, split = "_")[[1]]
  enh <- ar[1]
  gene <- ar[2]
  v <- jaccard_pseudotime(x = gene.mat[gene,], y = enh.mat[enh,], 
                          z = pt.cells, nBins = 5, div.method = "min", ifPlot = F)
  v[[1]]
}, mc.cores = detectCores()))
heatmap(t(scale(t(link.mat[1:100,]))), Colv = NA)
dim(link.mat)

spearman.mat <- pbmclapply(1:nrow(link.mat), function(i) {
  cor.test(link.mat[i, ], 1:ncol(link.mat), method = "spearman")$estimate
},mc.cores = detectCores()) %>% unlist
range(spearman.mat)
names(spearman.mat) <- seq_along(spearman.mat)
spearman.mat <- spearman.mat[!is.na(spearman.mat)] %>% sort
head(spearman.mat)



# Divide links into three groups
link.sets <- data.frame(cells = names(spearman.mat), bin = cut(spearman.mat, 3)) %>% 
         split(., .$bin) %>% sapply(., "[[", "cells") %>% sapply(., function(x) {
           u.links[as.numeric(x)]
         })
names(link.sets)
length(links)
link.fr <- Reduce("rbind", pbmclapply(links, function(x) {
  sapply(link.sets, function(y) {
    length(intersect(x, y)) / length(x)
  })
}, mc.cores = detectCores()))
rownames(link.fr) <- seq_along(links)
apply(link.fr, 2, range)

down.link.ids <- sort(link.fr[, 1], decreasing = T) %>% head(n = 5) %>% names %>% as.numeric
down.link.ids

mid.link.ids <- sort(link.fr[, 2], decreasing = T) %>% head(n = 5) %>% names %>% as.numeric
mid.link.ids

up.link.ids <- sort(link.fr[, 3], decreasing = T) %>% head(n = 5) %>% names %>% as.numeric
up.link.ids



# Build matrix for the down links
down.mat <- Reduce("rbind", lapply(down.link.ids, function(x) {
  xx <- intersect(links[[x]], link.sets[[1]]) %>% strsplit(., split = "_")
  Reduce("rbind",pbmclapply(xx, function(y) {
    enh <- y[1]
    gene <- y[2]
    v <- jaccard_pseudotime(x = gene.mat[gene,], y = enh.mat[enh,], 
                            z = pt.cells, nBins = 20, div.method = "min", ifPlot = F)
    v[[1]]
  }, mc.cores = detectCores())) %>% t %>% scale %>% t %>% apply(., 2, mean)
}))
rownames(down.mat) <- down.link.ids



# Build matrix for the mid links
mid.mat <- Reduce("rbind", lapply(mid.link.ids, function(x) {
  xx <- intersect(links[[x]], link.sets[[2]]) %>% strsplit(., split = "_")
  Reduce("rbind",pbmclapply(xx, function(y) {
    enh <- y[1]
    gene <- y[2]
    v <- jaccard_pseudotime(x = gene.mat[gene,], y = enh.mat[enh,], 
                            z = pt.cells, nBins = 20, div.method = "min", ifPlot = F)
    v[[1]]
  }, mc.cores = detectCores())) %>% t %>% scale %>% t %>% apply(., 2, mean)
}))
rownames(mid.mat) <- mid.link.ids



# Build matrix for the up links
up.mat <- Reduce("rbind", lapply(up.link.ids, function(x) {
  xx <- intersect(links[[x]], link.sets[[3]]) %>% strsplit(., split = "_")
  Reduce("rbind",pbmclapply(xx, function(y) {
    enh <- y[1]
    gene <- y[2]
    v <- jaccard_pseudotime(x = gene.mat[gene,], y = enh.mat[enh,], 
                            z = pt.cells, nBins = 20, div.method = "min", ifPlot = F)
    v[[1]]
  }, mc.cores = detectCores())) %>% t %>% scale %>% t %>% apply(., 2, mean)
}))
rownames(up.mat) <- up.link.ids



# Combine the matrices
ht.mat <- Reduce("rbind", list(down.mat, mid.mat, up.mat))
heatmap(ht.mat, Rowv = NA, Colv = NA)
range(ht.mat)
rownames(ht.mat) <- sapply(cts.en.regs, "[[", "TF")[as.numeric(rownames(ht.mat))] %>% 
  strsplit(., split = "::") %>% sapply(., function(x) {tail(x, n = 1)}) %>% 
  tolower %>% capitalize
p.heatmap <- get_simple_heatmap(mat = ht.mat[, -c(5, 8)], name = "Z score of regulatory strength", 
                   bottom = -2, interium = 0, up = 0.4, labels = c("-3.4", "1.2"),
                   low = "white", middle = "#abedf5", high = "#31E1F7")
save_image(p = p.heatmap, paste0(image.dir, "IN_heatmap.png"), 
           width = 2000, height = 2000)



# Pseudotime plot
p.pt <- plot_cells(cds,
                   color_cells_by = "pseudotime", label_roots = F,
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=1.5)
save_image(p = p.pt, paste0(image.dir, ct, "_UMAP_pseudotimes.png"), 
           width = 2500, height = 2000)


# Get enhancers relevant to down regulation
down.link.ids
down.enh <- Reduce("c", lapply(down.link.ids, function(x) {
  xx <- intersect(links[[x]], link.sets[[1]])
  strsplit(xx, split = "_") %>% sapply(., "[[", 1)
})) %>% table %>% sort(., decreasing = T)
head(down.enh)
down.id <- 1
down.en.regs <- cts.en.regs[down.link.ids][sapply(links[down.link.ids], function(x) {
  names(down.enh[down.id]) %in% (strsplit(x, split = "_") %>% sapply(., "[[", 1) %>% unique)
}) %>% unlist]
sapply(down.en.regs, "[[", "TF")
down.genes <- sapply(links[down.link.ids], function(x) {
  xx <- intersect(x, link.sets[[1]])
  (strsplit(xx, split = "_") %>% sapply(., "[[", 2))[grepl(names(down.enh[down.id]), xx)]
}) %>% Reduce("union", .)
down.genes



# Plotting for down
down.dot <- jaccard_pseudotime(x = gene.mat[down.genes[1],], y = enh.mat[names(down.enh)[down.id],], 
                               eq.x = 1.0, eq.y = 0.15, 
                   rsquare.x = 1.0, rsquare.y = 0.35,
                   z = pt.cells, nBins = 20, div.method = "min")
save_image(p = down.dot, path = paste0(image.dir, ct, "_down_dotplot.png"), width = 1000, 
           height = 1000)



# Get enhancers relevant to up regulation
up.link.ids
up.enh <- Reduce("c", lapply(up.link.ids, function(x) {
  xx <- intersect(links[[x]], link.sets[[3]])
  strsplit(xx, split = "_") %>% sapply(., "[[", 1)
})) %>% table %>% sort(., decreasing = T)
head(up.enh)
up.id <- 7
up.en.regs <- cts.en.regs[up.link.ids][sapply(links[up.link.ids], function(x) {
  names(up.enh[up.id]) %in% (strsplit(x, split = "_") %>% sapply(., "[[", 1) %>% unique)
}) %>% unlist]
sapply(up.en.regs, "[[", "TF")
up.genes <- sapply(links[up.link.ids], function(x) {
  xx <- intersect(x, link.sets[[3]])
  (strsplit(xx, split = "_") %>% sapply(., "[[", 2))[grepl(names(up.enh[up.id]), xx)]
}) %>% Reduce("union", .)
up.genes



# Plotting for up
up.dot <- jaccard_pseudotime(x = gene.mat[up.genes[1],], y = enh.mat[names(up.enh)[up.id],], 
                             eq.x = 1.0, eq.y = 0.15, 
                   rsquare.x = 1.0, rsquare.y = 0.35,
                   z = pt.cells, nBins = 20, div.method = "min")
save_image(p = up.dot, path = paste0(image.dir, ct, "_up_dotplot.png"), width = 1000, 
           height = 1000)



save.image(paste0(scratch.dir, "24_IN_analysis.RData"))
