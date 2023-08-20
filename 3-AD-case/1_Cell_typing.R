############################################
#                                          #  
# Build heatmaps using Shane's heatmaps    #
#                                          #  
############################################



# Parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_3_AD/Rdata/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_3_AD/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_3_AD/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
setwd(R.dir)



# Sources
source(paste0(tool.dir, "visual_tools.R"))
load("5_AD_cell_typing.RData")



# Load data
obj <- qs::qread(paste0(R.dir, "Obj_clustered.qsave"))


dim(obj)
unique(obj$Time)
dim(obj@assays$RNA@scale.data)
dim(obj@assays$Activity@scale.data)
rownames(obj@assays$Activity@scale.data) # only 40 genes
markers <- read_excel(paste0(table.dir, "Markers_from_Shane.xlsx"), sheet = 1)
dim(markers)
markers



# Convert upper cases into lower ones
markers <- apply(markers, 2, function(x) {
  capitalize(tolower(x))
}) %>% as.data.frame



# Build heatmap for gene expression
obj$celltypes <- factor(obj$celltypes, levels = 0:(length(unique(obj$celltypes)) - 1))
get_marker_heatmap(obj = obj, markers = markers, 
                   path = paste0(image.dir, "Heatmaps_expression_activity.png"), 
                   ct.label = "celltypes")




# Load annotation from Shane
annot.df <- read.csv(paste0(table.dir, "Cell_type_annot_Shane_exclude_unlabeled.csv"))
dim(annot.df)
dim(obj)
head(annot.df)
length(unique(annot.df[, 1]))
setdiff(unique(annot.df[, 1]), colnames(obj)) %>% length # 12085
setdiff(colnames(obj), unique(annot.df[, 1]))
colnames(annot.df) <- c("barcode", "cell.type")
bc.df <- annot.df[annot.df$barcode %in% colnames(obj) & 
                    annot.df$cell.type != "",]
dim(bc.df)
bc.ct <- bc.df$cell.type
names(bc.ct) <- bc.df$barcode
head(bc.ct)
obj <- AddMetaData(obj, metadata = bc.ct, col.name = "Shane.cell.type")
table(obj$Shane.cell.type)
unique(obj$Shane.cell.type)
qs::qsave(obj, paste0(R.dir, "Object_Shane_annot.qsave"))
p.umap <- get_simple_UMAP(obj, group.by = "Shane.cell.type")
p.umap
qs::qsave(p.umap, paste0(R.dir, "UMAP_Shane_annot.qsave"))
save_image(p = p.umap, path = paste0(image.dir, "UMAP_Shane_annot.png"), 
           width = 2500, height = 2500)
obj <- FindClusters(obj, algorithm = 3, graph.name = "wsnn", resolution = 0.2)
cell.type <- obj$seurat_clusters
levels(cell.type) <- c(
  "Oligo", "EX", "IN", "AG", "EX", "IN", # 5
  "EX", "MG", "OPC", "IN", "IN", "AG", "IN", "Endo & Pericyte", "EX", "IN", # 15
  "EX", "EX", "Endo & Pericyte", "EX", "EX", "IN", "EX", # 22
  "AG", "EX", "OPC", "EX"
)
obj <- AddMetaData(obj, col.name = "Shane.cell.type", metadata = cell.type)
obj$celltypes <- obj$seurat_clusters
get_marker_heatmap(obj = obj, markers = markers, 
                   path = paste0(image.dir, "Heatmaps_expression_activity_seurat_clusters.png"), 
                   ct.label = "celltypes")
obj$celltypes <- obj$Shane.cell.type
get_marker_heatmap(obj = obj, markers = markers, 
                   path = paste0(image.dir, "Heatmaps_expression_activity_Shane_annot.png"), 
                   ct.label = "celltypes")
sink(paste0(table.dir, "Celltype_clusters.txt"))
print(ct.cl)
sink()



save.image("5_AD_cell_typing.RData")
