###################################################
#                                                 #
# Select ChIP-seq data to validate the identified #
# TFs and TF-region relations                     #
#                                                 #
###################################################


# General libraries
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
  "igraph", 
  "sparseMatrix", 
  "readxl", 
  "writexl"
)
libraries(libs)


# Parameters
data.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/"


# Load data
message ("\nBegan reading the Excel sheets one by one ...")
sheet_names <- excel_sheets(paste0(code.dir, "ATAC-GEX-TF-occupancy-validation.xlsx"))
message (length(sheet_names), " sheets exist in this Excel file")
lost.files <- c(
  "K562-hg38-ChIP-seq", 
  "BJ-hg38-ChIP-seq", 
  "GM12878-hg19-ChIP-seq", 
  "A549-hg38-ChIP-seq"
)
selected.ll <- list()
for(i in seq_along(sheet_names)) {
  
  # Read the sheet, skipping the first row
  message ("----> Reading the sheet ", sheet_names[i], " ...")
  data <- read_excel(paste0(code.dir, "ATAC-GEX-TF-occupancy-validation.xlsx"), 
                     sheet = sheet_names[i], skip = 1) %>% filter (`File Format` == "bed")
  # if (grepl(".bed.gz", data$`Target label`[1])) {
  #   target.chipseq <- split(data, f = data$`Batch`)
  # } else {
  #   target.chipseq <- split(data, f = data$`Target label`)
  # }
  tf.vec <- sapply(1:nrow(data), function(j) {
    ifelse (grepl(".bed.gz", data$`Target label`[j]), return(data$`Batch`[j]), 
            return(data$`Target label`[j]))
  })
  data <- cbind(data, TF = tf.vec)
  target.chipseq <- split(data, f = data$TF)
  selected <- lapply(seq_along(target.chipseq), function(j) {
    message ("--------> Processing ChIP-seq data of target ", names(target.chipseq)[j], " ...")
    chipseq.ll <- target.chipseq[[j]]
    chipseq.ll[order(chipseq.ll$`File size`, decreasing = TRUE), ] %>% head(n = 1)
  }) %>% Reduce("rbind", .)
  selected.ll[[i]] <- selected
  message ("----> Finished reading the sheet ", sheet_names[i], "\n")
  
  
  # Download dataset
  # if (!sheet_names[i] %in% lost.files) {
  #   next
  # }
  # if (i != 9) {
  #   next
  # }
  message ("\n----> Downloading datasets listed in the sheet ", sheet_names[i], " ...")
  dir.create(paste0(data.dir, sheet_names[i]))
  lapply(1:nrow(selected), function(j) {
    message ("--------> Downloading ChIP-seq data of ID ", j, " ...")
    if (grepl("bed.gz", selected$`Download URL`[j])) {
      url <- selected$`Download URL`[j]
      target <- selected$`Target label`[j]
    } else {
      url <- selected$`Target label`[j]
      target <- selected$Batch[j]
    }
    bed.file <- strsplit(url, split = "/")
    suffix <- bed.file[[1]][length(bed.file[[1]])]
    messag <- try (download.file(url = paste0("https://www.encodeproject.org", url), 
                  destfile = paste0(data.dir, sheet_names[i], "/", 
                                    target, "_", suffix) ) )
    if (class(messag) == "try-error") {
      message ("--------> Failed to download file for ", target)
    }
    message ("--------> Finished downloading ChIP-seq data of target ", target, "\n")
  })
}
names(selected.ll) <- sheet_names
write_xlsx(selected.ll, paste0(code.dir, "Selected-ATAC-GEX-TF-occupancy-validation.xlsx"))
message ("Finished reading the Excel sheets one by one\n")


save.image(paste0(data.dir, "1-Select-TF-ChiPseq.RData"))
# load("/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/1-Select-TF-ChiPseq.RData")