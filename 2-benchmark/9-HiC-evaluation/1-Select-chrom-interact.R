###################################################
#                                                 #
# Select Hi-C data to validate the identified     #
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
  "Matrix", 
  "readxl", 
  "writexl"
)
libraries(libs)


# Parameters
data.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/Chrom-interaction-validation/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Evaluation/Chrom-interaction-validation/"


# Load data
message ("\nBegan reading the Excel sheets one by one ...")
sheet_names <- excel_sheets(paste0(code.dir, "ATAC-GEX-chrom-interaction-validation.xlsx"))
message (length(sheet_names), " sheets exist in this Excel file")
#lost.file <- "A549-hg38-HiC"
selected.ll <- list()
for(i in seq_along(sheet_names)) {
  
  # Read the sheet, skipping the first row
  message ("----> Reading the sheet ", sheet_names[i], " ...")
  data <- read_excel(paste0(code.dir, "ATAC-GEX-chrom-interaction-validation.xlsx"), 
                     sheet = sheet_names[i], skip = 1) %>% filter (`File Format` == "hic")
  # if ("CTCF" %in% data$`Target label`) {
  #   target.hic <- list(CTCF = split(data, f = data$`Target label`) %>% `[[` ("CTCF"))
  # } else if ("POLR2A" %in% data$`Target label`) {
  #   target.hic <- list(POLR2A = split(data, f = data$`Target label`) %>% `[[` ("POLR2A"))
  # } else {
  #   target.hic <- list("NA" = data)
  # }
  # if (!"POLR2A" %in% data$`Target label`) {
  #   next
  # }
  #target.hic <- list(POLR2A = split(data, f = data$`Target label`) %>% `[[` ("POLR2A"))
  target.hic <- list(data)
  selected <- lapply(seq_along(target.hic), function(j) {
    message ("--------> Processing Hi-C data of target ", names(target.hic)[j], " ...")
    hic.ll <- target.hic[[j]]
    hic.ll[order(hic.ll$`File size`, decreasing = TRUE), ] %>% head(n = 1)
  }) %>% Reduce("rbind", .)
  selected.ll[[i]] <- selected
  message ("----> Finished reading the sheet ", sheet_names[i], "\n")
  
  
  # Download dataset
  # if (sheet_names[i] != lost.file) {
  #   next
  # }
  message ("\n----> Downloading datasets listed in the sheet ", sheet_names[i], " ...")
  dir.create(paste0(data.dir, sheet_names[i]))
  lapply(1:nrow(selected), function(j) {
    message ("--------> Downloading Hi-C data of target ", selected$`Target label`[j], " ...")
    hic.file <- strsplit(selected$`Download URL`[j], split = "/")
    suffix <- hic.file[[1]][length(hic.file[[1]])]
    messag <- try(download.file(url = paste0("https://www.encodeproject.org", selected$`Download URL`[j]), 
                  destfile = paste0(data.dir, sheet_names[i], "/", 
                                    selected$`Target label`[j], "_", suffix) ) )
    if (class(messag) == "try-error") {
      message ("--------> Failed to download file for ", selected$`Target label`[j])
    }
    message ("--------> Finished downloading Hi-C data of target ", selected$`Target label`[j], "\n")
  })
}
names(selected.ll) <- sheet_names
write_xlsx(selected.ll, paste0(code.dir, "Selected-ATAC-GEX-chrom-interaction-validation.xlsx"))
message ("Finished reading the Excel sheets one by one\n")


save.image(paste0(data.dir, "1-Select-chrom-interact.RData"))
# load("/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/Chrom-interaction-validation/1-Select-chrom-interact.RData")