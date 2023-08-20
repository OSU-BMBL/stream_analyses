#####################################################
#                                                   #
# Downsample data to incorporate 2,000 cells,       #
# 10,000 genes, and the first half variable peaks   #
#                                                   #
#####################################################


# Task list
# 1. Remove the non-gene GEX features, which begin with "AC", "AL", or "LINC" 
#    as well as the ones on mitochondrion DNA sequences, i.e., the ones begin 
#    with "MT" or "mt"
# 2. Remove peaks on non-standard chromosomes
# 3. Select the top-ranked 10,000 variable genes
# 4. Choose the top-ranked half variable peaks
# 5. Randomly select 2,000 cells



# Libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(easypackages)
libs <- c(
  "qs", 
  "hdf5r",
  "Seurat",
  "Signac", 
  "pbmcapply", 
  "pbapply", 
  "parallel",
  "EnsDb.Mmusculus.v79", 
  "EnsDb.Hsapiens.v86",
  "dplyr", 
  "ggplot2", 
  "ggpubr", 
  "stringr"
)
libraries(libs)


# Parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Downsampling/"
parent.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/"
data.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/"
dir.create("/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/")
# org <- "hg38"
setwd(data.dir)
getwd()


# List folders composed of the whole data
dir.list <- setdiff(list.dirs(paste0(data.dir, "../"), full.names = FALSE, 
                              recursive = FALSE), "Downsampled-data")
length(dir.list)
head(dir.list)


# Downsampling each dataset
pblapply(dir.list, function(x) {
  message (x)
  dir.create(paste0(data.dir, x))
  rda.list <- list.files(paste0(parent.dir, x), pattern = ".qsave")
  frag.list <- list.files(paste0(parent.dir, x), pattern = "_fragments.tsv.gz")

  
  lapply(rda.list, function(y) {
    message ("----> ", y)
    org <- str_extract(y, regex("hg19|hg38|mm10|mm9", ignore_case = TRUE))
    ifelse(grepl("^mm", org), mt.prefix <- "mt-", mt.prefix <- "MT-")
    obj <- qs::qread(paste0(parent.dir, x, "/", y))
    message ("----> Cells: ", ncol(obj), "\n", 
             "----> Genes: ", nrow(obj[['RNA']]), "\n", 
             "----> Peaks: ", nrow(obj[['ATAC']]), "\n")
    
    
    # Randomly select cells
    obj.save <- obj
    n.cells <- min(1000, ncol(obj))
    obj <- subset(obj, cells = colnames(obj)[sample(x = 1:ncol(obj), size = n.cells)])
    message ("----> ", ncol(obj), " cells were selected after downsampling")
    
    
    # Filter out non-genes or MT-genes
    # genes <- setdiff(rownames(obj[['RNA']]), 
    #                  intersect(grep("^AC|^AL|^LINC", rownames(obj[['RNA']]), value = TRUE), 
    #                            grep("\\.", rownames(obj[['RNA']]), value = TRUE)) )
    genes <- setdiff(rownames(obj[['RNA']]), 
                     grep("(^AC|^AL|^LINC)[0-9]+", rownames(obj[['RNA']]), 
                          value = TRUE) )
    message ("----> ", length(genes), " genes after filtering out primer names")
    genes <- grep(mt.prefix, genes, invert = TRUE, value = TRUE)
    message ("----> ", length(genes), " genes after filtering out MT-genes")
    
    
    # Select the top-ranked variable genes
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst")
    rna.meta <- obj[['RNA']]@meta.features
    rna.meta <- rna.meta[order(rna.meta$vst.variance.standardized, decreasing = TRUE),]
    genes <- rna.meta[rownames(rna.meta) %in% genes,] %>% head(n = 10000) %>% rownames
    message ("----> ", length(genes), " genes after selecting highly variable genes")
    
    
    # Filter out peaks on non-standard chromosomes
    peaks <- rownames(obj[['ATAC']])
    peaks <- grep("chr([0-9]|X|Y)", peaks, value = TRUE)
    message ("----> ", length(peaks), " genes after filtering out peaks on non-standard chromosomes")
    
    
    # Select the top-ranked variable peaks
    DefaultAssay(obj) <- "ATAC"
    obj <- FindTopFeatures(obj)
    atac.meta <- obj[['ATAC']]@meta.features
    atac.meta <- atac.meta[order(atac.meta$percentile, decreasing = TRUE),]
    peaks <- atac.meta[rownames(atac.meta) %in% peaks,] %>% head(n = 10000) %>% rownames
    message ("----> ", length(peaks), " peaks after selecting highly variable peaks")
    
    
    # Create the subsetted Seurat object
    obj.sub <- subset(obj, features = c(genes, peaks))
    message ("----> ", nrow(obj.sub[['RNA']]), " genes, ", 
             nrow(obj.sub[['ATAC']]), " peaks, and ", ncol(obj.sub), 
             " cells were retained after downsampling")
    
    
    # Save the Seurat object and copy fragments file
    qs::qsave(obj.sub, paste0(data.dir, x, "/", y))
    if (length(frag.list) > 0) {
      system(paste0("cp ", data.dir, "../", x, "/", gsub(".qsave", "_fragments.tsv.gz ", y), 
                    data.dir, x))
    }
    message ("----> Finished copying ", gsub(".qsave", "_fragments.tsv.gz ", y, "\n\n"))
  })
})