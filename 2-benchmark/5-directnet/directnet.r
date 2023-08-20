library("tibble")
library(rlang)
library(dplyr)#version should be 1.0.7
library(presto)
library("JASPAR2016")
library("spatstat.sparse")
library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(ggplot2)
options(stringsAsFactors = FALSE)
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
library("qs")
library(BSgenome.Hsapiens.UCSC.hg38)

# we used default cutoff from Seurat website in the function "Findmarkers" logfc.

directnet_run<-function(file.path,output_path,output_name){
  #geneinfo
  genome.info <- read.table(file = "/fs/ess/PCON0022/guoqi/Yang/Stream/hg38.promoter.regions.txt")
  names(genome.info) <- c("Chrom","Starts","Ends","genes")
  genes <- lapply(genome.info$genes, function(x) strsplit(x,"[|]")[[1]][1])
  genes <- lapply(genes, function(x) strsplit(x,"[.]")[[1]][1])
  genes <- unlist(genes)
  genome.info$genes <- genes
  unik <- !duplicated(genes)# filter out different transcript
  genome.info <- genome.info[unik,]
  
  #load seurat object
  examp_3_3_3 <- qread(file.path)
  #process rna
  DefaultAssay(examp_3_3_3) <- "RNA"
  examp_3_3_3 <- SCTransform(examp_3_3_3, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  #process atac
  DefaultAssay(examp_3_3_3) <- "ATAC"
  examp_3_3_3 <- RunTFIDF(examp_3_3_3)
  examp_3_3_3 <- FindTopFeatures(examp_3_3_3, min.cutoff = 'q0')
  examp_3_3_3 <- RunSVD(examp_3_3_3)
  examp_3_3_3 <- RunUMAP(examp_3_3_3, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  #integrate
  examp_3_3_3 <- FindMultiModalNeighbors(examp_3_3_3, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  examp_3_3_3 <- RunUMAP(examp_3_3_3, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  examp_3_3_3 <- FindClusters(examp_3_3_3, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  #identify markers
  DefaultAssay(examp_3_3_3)<-"SCT"
  #qsave(examp_3_3_3,"./Directnet/examp_3_3.qs")
  #examp_3_3_3<-qread("./Directnet/examp_3_3.qs")
  markers_all <- presto:::wilcoxauc.Seurat(X = examp_3_3_3, group_by = 'wsnn_res.0.8', assay = 'data', seurat_assay = 'SCT')
  markers_all <- markers_all[which(markers_all$auc > 0.5), ,drop = FALSE ]
  markers_all <- markers_all[which(markers_all$padj < 0.05), ,drop = FALSE]
  markers <- data.frame(gene = markers_all$feature, group = markers_all$group)
  c <- unique(markers$group)
  marker_list <- list()
  for (i in 1:length(c)) {
    marker1<- markers_all[markers$group == c[i],]
    marker_list[[i]] <- as.character(marker1$feature[marker1$auc > 0.5])
  }
  markers_groups <- unique(unlist(marker_list))
  markers_groups <- lapply(markers_groups, function(x) strsplit(x,"[.]")[[1]][1])
  markers_groups <- unique(unlist(markers_groups))
  #network inference
  examp_3_3_3 <- Run_DIRECT_NET(examp_3_3_3, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap=0.5, 
                                size_factor_normalize = FALSE, genome.info = genome.info, focus_markers =markers_all$feature)
  direct.net_result <- Misc(examp_3_3_3, slot = 'direct.net')
  direct.net_result <- as.data.frame(do.call(cbind,direct.net_result)) # links for markers
  direct.net_result$function_type <- gsub("HF","HC",direct.net_result$function_type)
  direct.net_result$function_type <- gsub("Rest","MC",direct.net_result$function_type)
  direct.net_result$function_type <- gsub("LF","LC",direct.net_result$function_type)

  DefaultAssay(examp_3_3_3) <- 'ATAC'
  focused_markers <- markers
  groups <- unique(focused_markers$group)
  da_peaks_list <- list()
  Idents(examp_3_3_3)<-examp_3_3_3$wsnn_res.0.8
  for (i in 1:length(groups)) {
    print(i)
    da_peaks <- FindMarkers(
      object = examp_3_3_3,
      min.pct = 0.1,
      logfc.threshold = 0.25,
      ident.1 = groups[i],
      group.by = "wsnn_res.0.8",
      test.use = 'LR',
      only.pos = TRUE
    )
    da_peaks_list[[i]] <- da_peaks
  }
  # CRE-gene connections
  CREs_Gene <- generate_CRE_Gene_links(direct.net_result, markers = focused_markers)
  # Find focused CREs which is overlapped with DA
  Focused_CREs <- generate_CRE(L_G_record = CREs_Gene$distal, P_L_G_record = CREs_Gene$promoter, da_peaks_list)
  # detect TFs for distal CREs
  L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$distal, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = focused_markers)
  # detect TFs for Promoters
  P_L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$promoter, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = focused_markers)
  #Detect CRE-TF connections
  network_links <- generate_links_for_Cytoscape(L_G_record = Focused_CREs$L_G_record, L_TF_record, P_L_G_record = Focused_CREs$P_L_G_record, P_L_TF_record,groups)
  #Ouput Node attribute
  Node_attribute <- generate_node_for_Cytoscape(network_links,markers = focused_markers)
  list_output<-list()
  list_output[[1]]<-network_links
  list_output[[2]]<-Node_attribute
  names(list_output)<-c("link","node")
  qsave(list_output,paste(output_path,output_name,sep="/"))
  }


#------------------------------------------------------output

#load data (parameters)
args <- commandArgs(TRUE)
file.path <- args[1]
file.name <- strsplit(file.path, split = "\\/") %>% `[[` (1) %>% tail(n = 1)
data.dir <- paste0("/fs/ess/PCON0022/guoqi/Yang/Stream/Directnet/output/", 
                   gsub(".qsave", "", file.name))
#dir.create(data.dir)
#org <- strsplit(file.name, split = "_") %>% `[[` (1) %>% `[` (2)


# Run DIRECTNET
setwd(data.dir)
output_name=paste0(file.name,"_directnet_out.qs")
directnet_run(file.path = file.path,output_path = data.dir,output_name=output_name)

