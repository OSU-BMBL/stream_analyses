#####################################################
#                                                   #
#                 Infer trajectories                #
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
#     15. monocle3_TI : run monocle3 to infer trajectory      #
#                                                             #
###############################################################


# Input : 
# 1. red.dim : whether use the already calculated dimensions


monocle3_TI <- function(obj = NULL, assay = "RNA", red.dim = F, 
                        color_cells_by = "cluster",
                        time.label = NULL, use_partition = F, 
                        num_dim = 50, dir = "./", prefix = "tmp",
                        reduction_method = "UMAP", ciliated_genes = NULL, 
                        min.ncells = 10, min_expr = 0.1, alignment_group = NULL, 
                        residual_model_formula_str = NULL) {
  
  # Parameters
  message ("To infer trajectory for a transcriptomics dataset composed of ", 
           ncol(obj), " cells ...")
  
  
  # Load data into monocle
  require(monocle3)
  require(dplyr)
  message ("Coverting the ", assay, " assay of a Seurat object into monocle CellDataSet object ...\n", 
           "Using negative binomial distribution for read counts ...")
  gene_annotation <- cbind(rownames(obj[[assay]]@meta.features), 
                           obj[[assay]]@meta.features)
  colnames(gene_annotation)[1] <- "gene_short_name"
  gene_annotation <- cbind(gene_annotation, 
                           num_cells_expressed = apply(GetAssayData(obj, assay = assay, 
                                                                    slot = "counts") > 0, 
                                                       1, sum))
  cds <- new_cell_data_set(GetAssayData(obj, assay = assay, slot = "counts"), 
                           cell_metadata = obj@meta.data, 
                           gene_metadata = gene_annotation)
  
  
  # Pre-process the data and dimension reduction
  message ("Preprocessing the CellDataSet ...")
  if (!red.dim) {
    message ('Re-calculating dimension reduction ...')
    cds <- preprocess_cds(cds, num_dim = num_dim)
    if (!is.null(alignment_group)) {
      message ("Removing batch effects ...")
      cds <- align_cds(cds, alignment_group = alignment_group, 
                       residual_model_formula_str = residual_model_formula_str)
    }
    message ("Performing dimensionality reduction ...")
    cds <- reduce_dimension(cds, reduction_method = reduction_method)
  } else {
    message ('Using the dimension reduction results saved in the Seurat object ...')
    cds@int_colData@listData$reducedDims@listData$PCA <- obj[['pca']]@cell.embeddings
    cds@int_colData@listData$reducedDims@listData$UMAP <- obj[['umap']]@cell.embeddings
  }
  
  
  # Visualize the results
  save_image(p = plot_cells(cds, label_groups_by_cluster = FALSE, 
                            color_cells_by = color_cells_by, 
                            reduction_method = reduction_method), 
             path = paste(dir, prefix, "_monocle3_UMAP.png"))
  
  
  # Visualize how individual genes vary along the trajectory
  if (!is.null(ciliated_genes)) {
    message ("Visualizing expression for genes: ")
    ciliated_genes
    save_image (p = plot_cells(cds,
                               genes = ciliated_genes,
                               label_cell_groups = FALSE,
                               show_trajectory_graph = FALSE), 
                path = paste0(dir, prefix, "_monocle3_expr.png"))
  }
  
  
  # Cluster your cells
  message ("Clustering cells ...")
  cds <- cluster_cells(cds)
  # if (!red.dim) {
  #   message ("Clustering cells using monocle3 ...")
  #   cds <- cluster_cells(cds)
  # } else {
  #   message ("Using the already calculated cell clusters ...")
  # }
  save_image(p = plot_cells(cds, color_cells_by = "partition"), 
             path = paste0(dir, prefix, "_monocle3_partitions.png"))
  
  
  # Learn the trajectory graph
  message ("Learning trajectory graph ...")
  cds <- learn_graph(cds, use_partition = use_partition)
  save_image(p = plot_cells(cds,
                            color_cells_by = color_cells_by,
                            label_groups_by_cluster = FALSE,
                            label_leaves = FALSE,
                            label_branch_points=FALSE), 
             path = paste0(dir, prefix, "_monocle3_trajectory_graph.png"))
  
  
  # Order the cells in pseudotime
  if (is.null(time.label)) {
    time.label <- color_cells_by 
  }
  message ("Ordering cells ...")
  save_image(p = plot_cells(cds,
                            color_cells_by = color_cells_by,
                            label_cell_groups=FALSE,
                            label_leaves=TRUE,
                            label_branch_points=TRUE,
                            graph_label_size=1.5), 
             path = paste0(dir, prefix, "_monocle3_colored_cells.png"))
  message ("Please choose root nodes manually via 'cds <- order_cells(cds)' or programmatically ...")
  
  
  cds
}


# Order cells alongside trajectory
monocle3_order_cells <- function(cds = NULL, root.method = "manual", 
                                 dir = "./", prefix = "tmp", time_bin = NULL, partition = F,
                                 time.label = "Time") {
  
  # Load data
  message ("Loading CellDataSet to order cells using ", root.method, " method to determine root ...")
  cds
  
  
  # Plotting for time
  message ("Generating a UMAP plot with time labels ...")
  save_image(p = plot_cells(cds,
                            color_cells_by = time.label,
                            label_cell_groups=FALSE,
                            label_leaves=TRUE,
                            label_branch_points=TRUE,
                            graph_label_size=1.5), 
             path = paste0(dir, prefix, "_times.png"))
  
  
  # Order cells
  if (root.method == "manual") {
    message ("Orderringf cells alongside trajectory manually ...")
    cds <- order_cells(cds)
  } else {
    # a helper function to identify the root principal points:
    message ("Ordering cells alongside trajectory programatically ...")
    cds <- order_cells(cds, root_pr_nodes = 
                         get_earliest_principal_node(cds, time_bin = time_bin, 
                                                     partition = partition, 
                                                     time.label = time.label))
  }
  
  
  # Plotting
  save_image(p = plot_cells(cds,
                            color_cells_by = "pseudotime",
                            label_cell_groups=FALSE,
                            label_leaves=FALSE,
                            label_branch_points=FALSE,
                            graph_label_size=1.5), 
             path = paste0(dir, prefix, "_monocle3_pseudotime.png"))
  message ("Returning the CellDataSet ...")
  
  
  cds
}


###############################################################
#                                                             #
# 1. run_hyper_test : run hyper geometric tests between two   #
#    lists of sets#                                           #
#                                                             #
###############################################################


# Input :
# x : the first list of sets
# y : the second list of sets
# n : the total number of all elements


run_hyper_test <- function(x = NULL, y = NULL, n = NULL) {
  
  # Check inputs
  message ("Comparing a list of ", length(x), " sets against a list of ", 
           length(y), " sets.\n")
  
  
  # Comparison
  library(pbmcapply)
  library(data.table)
  hyper.df <- as.data.frame(do.call("rbind", pbmclapply(seq_along(x), function(i) {
    xx <- x[[i]]
    Reduce("rbind", lapply(seq_along(y), function(j) {
      yy <- y[[j]]
      overlap <- length(intersect(xx, yy))
      c(i, j, overlap, length(xx), length(yy),
        phyper(overlap - 1, length(yy), n - length(yy), length(xx), 
               lower.tail = F))
    }))
  }, mc.cores = detectCores())))
  colnames(hyper.df) <- c("x", "y", "Overlap", "Query", "Target", "Pval")
  rownames(hyper.df) <- NULL
  message ("Identified ", nrow(hyper.df[hyper.df$Pval * length(y) < 0.05,]), 
           " significantly enriched query-target pairs out of ", 
           nrow(hyper.df), " after Bonferroni correction.\n")
  
  
  hyper.df
}


# Load data
obj <- qs::qread(paste0(R.dir, "Object_Shane_annot.qsave"))
Idents(obj) <- obj$month
table(obj$month)
obj.ls <- SplitObject(obj, split.by = "month")
qs::qsave(obj.ls, paste0(R.dir, "obj_list.qsave"))
sapply(obj.ls, ncol)
obj.2.5 <- obj.ls$`2.5`
obj.5.7 <- obj.ls$`5.7`
obj.13 <- obj.ls$`13`
en.reg.2.5 <- qs::qread(paste0(R.dir, "Extended_eGRNs_2.5.qsave"))
en.reg.5.7 <- qs::qread(paste0(R.dir, "Extended_eGRNs_5.7.qsave"))
en.reg.13 <- qs::qread(paste0(R.dir, "Extended_eGRNs_13.qsave"))



# 2.5 months
cts.res.2.5 <- run_hyper_test(x = sapply(en.reg.2.5, "[[", "cells"), 
                              y = sapply(split(obj.2.5@meta.data, f = obj.2.5@meta.data$celltypes), 
                                         rownames), n = ncol(obj.2.5))
dim(cts.res.2.5)
sig.cts.2.5 <- cts.res.2.5[cts.res.2.5$Pval < 0.05,]
sig.cts.2.5$y <- levels(obj.2.5$celltypes)[sig.cts.2.5$y]
dim(sig.cts.2.5)
head(sig.cts.2.5)
length(unique(sig.cts.2.5$x))
qs::qsave(sig.cts.2.5, paste0(R.dir, "sig_cts_2.5.qsave"))



# 5.7 months
cts.res.5.7 <- run_hyper_test(x = sapply(en.reg.5.7, "[[", "cells"), 
                              y = sapply(split(obj.5.7@meta.data, f = obj.5.7@meta.data$celltypes), 
                                         rownames), n = ncol(obj.5.7))
dim(cts.res.5.7)
sig.cts.5.7 <- cts.res.5.7[cts.res.5.7$Pval < 0.05,]
sig.cts.5.7$y <- levels(obj.5.7$celltypes)[sig.cts.5.7$y]
dim(sig.cts.5.7)
head(sig.cts.5.7)
length(unique(sig.cts.5.7$x))
qs::qsave(sig.cts.5.7, paste0(R.dir, "sig_cts_5.7.qsave"))



# 13 months
cts.res.13 <- run_hyper_test(x = sapply(en.reg.13, "[[", "cells"), 
                              y = sapply(split(obj.13@meta.data, f = obj.13@meta.data$celltypes), 
                                         rownames), n = ncol(obj.13))
dim(cts.res.13)
sig.cts.13 <- cts.res.13[cts.res.13$Pval < 0.05,]
sig.cts.13$y <- levels(obj.13$celltypes)[sig.cts.13$y]
dim(sig.cts.13)
head(sig.cts.13)
length(unique(sig.cts.13$x))
qs::qsave(sig.cts.13, paste0(R.dir, "sig_cts_13.qsave"))




# Infer trajectory for EX
Idents(obj) <- obj$celltypes
ex.obj <- subset(obj, ident = "EX")
dim(obj)
dim(ex.obj)
ex.cds <- monocle3_TI(obj = ex.obj, dir = image.dir, prefix = "ex", 
                       color_cells_by = "Time")
ex.cds <- monocle3_order_cells(cds = ex.cds, root.method = "programmatical", 
                                dir = image.dir, prefix = "ex", time_bin = "2.5 months", 
                                time.label = "Time", partition = T)
qs::qsave(ex.cds, paste0(R.dir, "ex_cds.qsave"))



# Infer trajectory for IN
in.obj <- subset(obj, ident = "IN")
dim(obj)
dim(in.obj)
in.cds <- monocle3_TI(obj = in.obj, dir = image.dir, prefix = "in", 
                      color_cells_by = "Time")
in.cds <- monocle3_order_cells(cds = in.cds, root.method = "programmatical", 
                               dir = image.dir, prefix = "in", time_bin = "2.5 months", 
                               time.label = "Time", partition = T)
qs::qsave(in.cds, paste0(R.dir, "in_cds.qsave"))



# Infer trajectory for Oligo
oligo.obj <- subset(obj, ident = "Oligo")
dim(obj)
dim(oligo.obj)
oligo.cds <- monocle3_TI(obj = oligo.obj, dir = image.dir, prefix = "oligo", 
                      color_cells_by = "Time")
oligo.cds <- monocle3_order_cells(cds = oligo.cds, root.method = "programmatical", 
                               dir = image.dir, prefix = "oligo", time_bin = "2.5 months", 
                               time.label = "Time", partition = T)
qs::qsave(oligo.cds, paste0(R.dir, "oligo_cds.qsave"))



# Infer trajectory for Astrocytes
ag.obj <- subset(obj, ident = "AG")
dim(obj)
dim(ag.obj)
ag.cds <- monocle3_TI(obj = ag.obj, dir = image.dir, prefix = "ag", 
                         color_cells_by = "Time")
ag.cds <- monocle3_order_cells(cds = ag.cds, root.method = "programmatical", 
                                  dir = image.dir, prefix = "ag", time_bin = "2.5 months", 
                                  time.label = "Time", partition = T)
qs::qsave(ag.cds, paste0(R.dir, "ag_cds.qsave"))



# Infer trajectory for Microglia
mg.obj <- subset(obj, ident = "MG")
dim(obj)
dim(mg.obj)
mg.cds <- monocle3_TI(obj = mg.obj, dir = image.dir, prefix = "mg", 
                      color_cells_by = "Time")
mg.cds <- monocle3_order_cells(cds = mg.cds, root.method = "programmatical", 
                               dir = image.dir, prefix = "mg", time_bin = "2.5 months", 
                               time.label = "Time", partition = T)
qs::qsave(mg.cds, paste0(R.dir, "mg_cds.qsave"))



# Infer trajectory for OPC
opc.obj <- subset(obj, ident = "OPC")
dim(obj)
dim(opc.obj)
opc.cds <- monocle3_TI(obj = opc.obj, dir = image.dir, prefix = "opc", 
                      color_cells_by = "Time")
opc.cds <- monocle3_order_cells(cds = opc.cds, root.method = "programmatical", 
                               dir = image.dir, prefix = "opc", time_bin = "2.5 months", 
                               time.label = "Time", partition = T)
qs::qsave(opc.cds, paste0(R.dir, "opc_cds.qsave"))



save.image(paste0(scratch.dir, "22_trajectory_celltypes.RData"))
