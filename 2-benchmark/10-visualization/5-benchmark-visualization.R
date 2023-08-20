################################################################
#                                                              #
#            Generate figures for benchmarking                 #
#                                                              #
################################################################


# General libraries
dyn.load(x = "/users/PAS1475/liyang/libs/hdf5_1.10.6/lib/libhdf5_hl.so.100")
library(stream2)
message ("Loaded libraries: ", paste(libs, collapse = ", "))


# Set parameters
code.dir <- "/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/"
parent.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/"
chipseq.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/TF-occupancy-validation/"
hic.dir <- "/fs/scratch/PCON0022/liyang/STREAM-revision/Evaluation/Chrom-interaction-validation/"
setwd(parent.dir)
getwd()


# Panels:
# 1. Chart of benchmarking pipeline
# 2. Numbers of TFs, regions, and genes per regulon for six methods
# 3. Bar plots of f-scores of TF coverage
# 4. Box plots of precision of region-to-gene linkages with P-values
# 5. Bar plots of f-scores of region-to-gene linkages
# 6. TF-target validation via TF perturbation


################################################################
#                                                              #
#            Obtain the representative data sets               #
#                                                              #
################################################################


benchTF.scores <- qs::qread(paste0(parent.dir, "TF-coverage-benchmark.qsave"))
benchTF.scores <- benchTF.scores[sapply(benchTF.scores, sum) > 0]
length(benchTF.scores)
benchTFRegion.scores <- qs::qread(paste0(parent.dir, "TF-to-region-benchmark.qsave"))
benchTFRegion.scores <- benchTFRegion.scores[, apply(benchTFRegion.scores, 2, sum) > 0]
benchTFRegion.all <- qs::qread(paste0(parent.dir, "TF-to-region-benchmark-all.qsave"))
ncol(benchTFRegion.scores)
benchRegionGene.scores <- qs::qread(paste0(parent.dir, "Region-to-gene-benchmark.qsave"))
benchRegionGene.scores <- benchRegionGene.scores[sapply(benchRegionGene.scores, sum) > 0]
length(benchRegionGene.scores)


bench.data <- Reduce("intersect", 
                       list(
                         names(benchTF.scores), 
                         colnames(benchTFRegion.scores), 
                         names(benchRegionGene.scores)
                       ))
bench.data


# TFs
benchTF.sel <- benchTF.scores[bench.data][sapply(benchTF.scores[bench.data], function(x) {
  which.max(x[, 3]) == 1
})]
length(benchTF.sel)


# TF-to-region
benchTfRegion.sel <- benchTFRegion.scores[, bench.data][, apply(benchTFRegion.scores[, bench.data], 2, function(x) {
  which.max(x) == 4
})]
ncol(benchTfRegion.sel)

benchTFRegion.final <- benchTFRegion.all[bench.data]
mean.fscore <- sapply(seq_along(benchTFRegion.final), function(i) {
  message (i)
  x <- benchTFRegion.final[[i]]
  split(x, f = x$method) %>% sapply(., function(xx) {
    mean(as.numeric(xx$f.score) )
  })
})
colnames(mean.fscore) <- names(benchTFRegion.final)


# Region-to-gene
benchRegionGene.sel <- benchRegionGene.scores[bench.data][sapply(benchRegionGene.scores[bench.data], function(x) {
  which.max(x[, 3]) == 1
})]
length(benchRegionGene.sel)


# Data for visualization
vis.data <- Reduce("intersect", 
       list(
         names(benchTF.sel), 
         colnames(benchTfRegion.sel), 
         names(benchRegionGene.sel)
       ))
vis.data
benchTF.sel[vis.data]
benchTfRegion.sel[, vis.data]
benchRegionGene.sel[vis.data]


################################################################
#                                                              #
#                    Load eRegulonbs/regulons                  #
#                                                              #
################################################################


#
# STREAM
# 
stream.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir, "/stream/", vis.data[i], "/eRegulons.qsave"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  qs::qread(paste0(parent.dir, "/stream/", vis.data[i], "/eRegulons.qsave"))
})
names(stream.regs) <- vis.data
qs::qsave(stream.regs, paste0(parent.dir, "stream_vis_regulons.qsave"))
message ("Number of datasets: ", length(stream.regs))


#
# SCENIC+
# 
scplus.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir, "/scenicplus/", vis.data[i], "/eRegulons.csv"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  read.csv(paste0(parent.dir, "/scenicplus/", vis.data[i], "/eRegulons.csv"), sep = ",")
})
names(scplus.regs) <- vis.data
length(scplus.regs)
head(scplus.regs[[5]])
sapply(scplus.regs, dim)
qs::qsave(scplus.regs, paste0(parent.dir, "scplus_vis_regulons.qsave"))
message ("Number of datasets: ", length(scplus.regs))


#
# SCENIC
#
scenic.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir,  "scenic/", vis.data[i], "/draft_grn.csv"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  read.csv(paste0(parent.dir, "scenic/", vis.data[i], "/draft_grn.csv"), sep = ",")
})
names(scenic.regs) <- vis.data
length(scenic.regs)
head(scenic.regs[[1]])
sapply(scenic.regs, dim)
qs::qsave(scenic.regs, paste0(parent.dir, "scenic_vis_regulons.qsave"))
message ("Number of datasets: ", length(scenic.regs))


#
# GLUE
#
scglue.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir, "/scglue/", vis.data[i], "/draft_grn.csv"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  grn <- read.csv(paste0(parent.dir, "/scglue/", vis.data[i], "/draft_grn.csv"), sep = ",")
  colnames(grn) <- c("TF", "Gene", "Strength")
  gene2peak <- read.csv(paste0(parent.dir, "/scglue/", vis.data[i], "/gene2peak.csv"), sep = ",")
  colnames(gene2peak) <- c("Gene", "Region", "Weight")
  common.genes <- intersect(unique(grn$Gene), unique(gene2peak$Gene) )
  merge(grn[grn$Gene %in% common.genes,, drop = FALSE], 
        gene2peak[gene2peak$Gene %in% common.genes,, drop = FALSE], by = "Gene")
})
names(scglue.regs) <- vis.data
length(scglue.regs)
head(scglue.regs[[1]])
sapply(scglue.regs, dim)
qs::qsave(scglue.regs, paste0(parent.dir, "scglue_vis_regulons.qsave"))
message ("Number of datasets: ", length(scglue.regs))


#
# DIRECT-NET
#
dinet.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir, "/directnet/", vis.data[i], "/", vis.data[i], 
                          ".qsave_directnet_out.qs"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  qs::qread(paste0(parent.dir, "/directnet/", vis.data[i], "/", vis.data[i], 
                         ".qsave_directnet_out.qs"))$link
})
names(dinet.regs) <- vis.data
length(dinet.regs)
head(dinet.regs[[1]])
sapply(dinet.regs, dim)
qs::qsave(dinet.regs, paste0(parent.dir, "directnet_vis_regulons.qsave"))
message ("Number of datasets: ", length(dinet.regs))


#
# Pando
#
pando.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir, "/pando/", vis.data[i], ".qs"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  dt <- qs::qread(paste0(parent.dir, "/pando/", vis.data[i], ".qs")) %>% NetworkModules
  dt@meta
})
names(pando.regs) <- vis.data
length(pando.regs)
head(pando.regs[[1]])
sapply(pando.regs, dim)
qs::qsave(pando.regs, paste0(parent.dir, "pando_vis_regulons.qsave"))
message ("Number of datasets: ", length(pando.regs))


#
# scMEGA
#
scmega.regs <- pblapply(seq_along(vis.data), function(i) {
  message ("----> Conversion for the ", i, "-th dataset: ", vis.data[i])
  if (!file.exists(paste0(parent.dir, "/scMEGA/", vis.data[i], ".csv"))) {
    message ("----> No eRegulon file on dataset: ", vis.data[i])
    return(NULL)
  }
  read.csv(paste0(parent.dir, "/scMEGA/", vis.data[i], ".csv"), sep = ",")
})
names(scmega.regs) <- vis.data
length(scmega.regs)
head(scmega.regs[[1]])
sapply(scmega.regs, dim)
qs::qsave(scmega.regs, paste0(parent.dir, "scmega_vis_regulons.qsave"))
message ("Number of datasets: ", length(scmega.regs))


################################################################
#                                                              #
#                       Numbers of TFs                         #
#                                                              #
################################################################


# Count number of TFs

# STREAM
nTfs.stream <- sapply(vis.data, function(x) {
  sapply(stream.regs[[x]], "[[", "TF") %>% unique %>% length
})
dfTfs.stream <- data.frame(
  No.TFs = nTfs.stream, 
  Method = rep("STREAM", length(nTfs.stream) )
)
dfTfs.stream


# SCENIC+
nTfs.scplus <- sapply(vis.data, function(x) {
  unique(scplus.regs[[x]][, "TF"]) %>% length
})
dfTfs.scplus <- data.frame(
  No.TFs = nTfs.scplus, 
  Method = rep("SCENIC+", length(nTfs.scplus) )
)
dfTfs.scplus


# SCENIC
nTfs.scenic <- sapply(vis.data, function(x) {
  unique(scenic.regs[[x]][, "TF"]) %>% length
})
dfTfs.scenic <- data.frame(
  No.TFs = nTfs.scenic, 
  Method = rep("SCENIC", length(nTfs.scenic) )
)
dfTfs.scenic


# GLUE
nTfs.scglue <- sapply(vis.data, function(x) {
  unique(scglue.regs[[x]][, "TF"]) %>% length
})
dfTfs.scglue <- data.frame(
  No.TFs = nTfs.scglue, 
  Method = rep("GLUE", length(nTfs.scglue) )
)
dfTfs.scglue


# DIRECT-NET
nTfs.dinet <- sapply(vis.data, function(x) {
  unique(unlist(dinet.regs[[x]][, "TF"]) ) %>% length
})
dfTfs.dinet <- data.frame(
  No.TFs = nTfs.dinet, 
  Method = rep("DIRECT-NET", length(nTfs.dinet) )
)
dfTfs.dinet


# Pando
nTfs.pando <- sapply(vis.data, function(x) {
  unique(unlist(pando.regs[[x]][, "tf"]) ) %>% length
})
dfTfs.pando <- data.frame(
  No.TFs = nTfs.pando, 
  Method = rep("Pando", length(nTfs.pando) )
)
dfTfs.pando


# scMEGA
nTfs.scmega <- sapply(vis.data, function(x) {
  unique(unlist(scmega.regs[[x]][, "tf"]) ) %>% length
})
dfTfs.scmega <- data.frame(
  No.TFs = nTfs.scmega, 
  Method = rep("scMEGA", length(nTfs.scmega) )
)
dfTfs.scmega


# Set colors
col.ls <- setNames(
  c("#C70039", "#614BC3", "#EA1179", "#1D5D9B", "#1A5D1A", "#FD8D14", "#068FFF"),
  c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA")
)


# Draw box plots to showcase No. of TFs
nTfs.df <- rbind(
  dfTfs.stream, 
  dfTfs.scplus, 
  dfTfs.scenic, 
  dfTfs.scglue, 
  dfTfs.dinet, 
  dfTfs.pando, 
  dfTfs.scmega
)
nTfs.df$Method <- factor(nTfs.df$Method, 
                         levels = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"))
library(ggplot2)
p.nTfs <- ggplot(nTfs.df, aes(x = Method, y = No.TFs, fill = Method)) + 
  geom_boxplot(lwd = 0.3) + scale_fill_manual(values = col.ls) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) + 
  theme(
    #axis.title = element_text(size = 24, color = "black"),        # Font settings for axis titles
    axis.text = element_text(size = 12, color = "black"),       # Font settings for axis text
    plot.title = element_text(size = 16, color = "black"), # Font settings for the plot title
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    panel.background = element_rect(fill = "white", linewidth = 0.3),
    legend.position = "none"
  )
nTfs.median <- split(nTfs.df, f = nTfs.df$Method) %>% sapply(., "[[", "No.TFs") %>% sapply(., median) %>% round
qs::qsave(p.nTfs, paste0(parent.dir, "boxplots_no_TFs.qsave"))
png(filename = paste0(parent.dir, "boxplots_no_TFs.png"), width = 300, height = 300, res = 150)
p.nTfs
dev.off()


################################################################
#                                                              #
#               Numbers of genes per regulon                   #
#                                                              #
################################################################


# STREAM
genesPerReg.stream <- sapply(vis.data, function(x) {
  sapply(stream.regs[[x]], "[[", "genes") %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.stream <- data.frame(
  No.genes.per.regulon = genesPerReg.stream, 
  Method = rep("STREAM", length(genesPerReg.stream) )
)
rownames(dfGenesPerReg.stream) <- NULL
dfGenesPerReg.stream


# SCENIC+
genesPerReg.scplus <- sapply(vis.data, function(x) {
  message (x)
  if (is.null(scplus.regs[[x]])) {
    return(NULL)
  }
  split(scplus.regs[[x]], f = scplus.regs[[x]]$TF) %>% sapply(., "[[", "Gene") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.scplus <- data.frame(
  No.genes.per.regulon = genesPerReg.scplus, 
  Method = rep("SCENIC+", length(genesPerReg.scplus) )
)
rownames(dfGenesPerReg.scplus) <- NULL
dfGenesPerReg.scplus


# SCENIC
genesPerReg.scenic <- sapply(vis.data, function(x) {
  if (is.null(scenic.regs[[x]])) {
    return(NULL)
  }
  split(scenic.regs[[x]], f = scenic.regs[[x]]$TF) %>% sapply(., "[[", "target") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.scenic <- data.frame(
  No.genes.per.regulon = genesPerReg.scenic, 
  Method = rep("SCENIC", length(genesPerReg.scenic) )
)
rownames(dfGenesPerReg.scenic) <- NULL
dfGenesPerReg.scenic


# GLUE
genesPerReg.scglue <- sapply(vis.data, function(x) {
  if (is.null(scglue.regs[[x]])) {
    return(NULL)
  }
  split(scglue.regs[[x]], f = scglue.regs[[x]]$TF) %>% sapply(., "[[", "Gene") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.scglue <- data.frame(
  No.genes.per.regulon = genesPerReg.scglue, 
  Method = rep("GLUE", length(genesPerReg.scglue) )
)
rownames(dfGenesPerReg.scglue) <- NULL
dfGenesPerReg.scglue


# DIRECT-NET
genesPerReg.dinet <- sapply(vis.data, function(x) {
  if (is.null(dinet.regs[[x]])) {
    return(NULL)
  }
  split(dinet.regs[[x]], f = dinet.regs[[x]]$TF) %>% sapply(., "[[", "Gene") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.dinet <- data.frame(
  No.genes.per.regulon = genesPerReg.dinet, 
  Method = rep("DIRECT-NET", length(genesPerReg.dinet) )
)
rownames(dfGenesPerReg.dinet) <- NULL
dfGenesPerReg.dinet


# Pando
genesPerReg.pando <- sapply(vis.data, function(x) {
  if (is.null(pando.regs[[x]])) {
    return(NULL)
  }
  split(pando.regs[[x]], f = pando.regs[[x]]$tf) %>% sapply(., "[[", "target") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.pando <- data.frame(
  No.genes.per.regulon = genesPerReg.pando, 
  Method = rep("Pando", length(genesPerReg.pando) )
)
rownames(dfGenesPerReg.pando) <- NULL
dfGenesPerReg.pando


# scMEGA
genesPerReg.scmega <- sapply(vis.data, function(x) {
  if (is.null(scmega.regs[[x]])) {
    return(NULL)
  }
  split(scmega.regs[[x]], f = scmega.regs[[x]]$tf) %>% lapply(., "[[", "gene") %>% lapply(., unique) %>% sapply(., length)
}) %>% unlist
dfGenesPerReg.scmega <- data.frame(
  No.genes.per.regulon = genesPerReg.scmega, 
  Method = rep("scMEGA", length(genesPerReg.scmega) )
)
rownames(dfGenesPerReg.scmega) <- NULL
dfGenesPerReg.scmega


# Draw box plots to showcase No. of genes per regulon
nGenesPerReg.df <- rbind(
  dfGenesPerReg.stream, 
  dfGenesPerReg.scplus, 
  dfGenesPerReg.scenic, 
  dfGenesPerReg.scglue, 
  dfGenesPerReg.dinet, 
  dfGenesPerReg.pando, 
  dfGenesPerReg.scmega
)
nGenesPerReg.df$Method <- factor(nGenesPerReg.df$Method, 
                         levels = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"))
p.nGenesPerReg.df <- ggplot(nGenesPerReg.df, aes(x = Method, y = No.genes.per.regulon, fill = Method)) + 
  geom_boxplot(lwd = 0.1) + scale_fill_manual(values = col.ls) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) + 
  theme(
    #axis.title = element_text(size = 24, color = "black"),        # Font settings for axis titles
    axis.text = element_text(size = 12, color = "black"),       # Font settings for axis text
    plot.title = element_text(size = 16, color = "black"), # Font settings for the plot title
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    panel.background = element_rect(fill = "white", linewidth = 0.3),
    legend.position = "none"
  )
nGenesPerReg.median <- split(nGenesPerReg.df, f = nGenesPerReg.df$No.genes.per.regulon) %>% 
  sapply(., "[[", "No.genes.per.regulon") %>% sapply(., median) %>% round
qs::qsave(p.nGenesPerReg.df, paste0(parent.dir, "boxplots_no_genes_per_regulon.qsave"))
png(filename = paste0(parent.dir, "boxplots_no_genes_per_regulon.png"), width = 300, height = 300, res = 150)
p.nGenesPerReg.df
dev.off()


################################################################
#                                                              #
#               Numbers of regions per regulon                 #
#                                                              #
################################################################


# STREAM
regionsPerReg.stream <- sapply(vis.data, function(x) {
  sapply(stream.regs[[x]], "[[", "peaks") %>% sapply(., length)
}) %>% unlist
dfRegionsPerReg.stream <- data.frame(
  No.regions.per.regulon = regionsPerReg.stream, 
  Method = rep("STREAM", length(regionsPerReg.stream) )
)
rownames(dfRegionsPerReg.stream) <- NULL
dfRegionsPerReg.stream


# SCENIC+
regionsPerReg.scplus <- sapply(vis.data, function(x) {
  message (x)
  if (is.null(scplus.regs[[x]])) {
    return(NULL)
  }
  split(scplus.regs[[x]], f = scplus.regs[[x]]$TF) %>% sapply(., "[[", "Region") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfRegionsPerReg.scplus <- data.frame(
  No.regions.per.regulon = regionsPerReg.scplus, 
  Method = rep("SCENIC+", length(regionsPerReg.scplus) )
)
rownames(dfRegionsPerReg.scplus) <- NULL
dfRegionsPerReg.scplus


# SCENIC
dfRegionsPerReg.scenic <- data.frame(
  No.regions.per.regulon = NA,
  Method = "SCENIC"
)
dfRegionsPerReg.scenic


# GLUE
regionsPerReg.scglue <- sapply(vis.data, function(x) {
  if (is.null(scglue.regs[[x]])) {
    return(NULL)
  }
  split(scglue.regs[[x]], f = scglue.regs[[x]]$TF) %>% sapply(., "[[", "Region") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfRegionsPerReg.scglue <- data.frame(
  No.regions.per.regulon = regionsPerReg.scglue, 
  Method = rep("GLUE", length(regionsPerReg.scglue) )
)
rownames(dfRegionsPerReg.scglue) <- NULL
dfRegionsPerReg.scglue


# DIRECT-NET
dfRegionsPerReg.dinet <- data.frame(
  No.regions.per.regulon = NA,
  Method = "DIRECT-NET"
)


# Pando
regionsPerReg.pando <- sapply(vis.data, function(x) {
  if (is.null(pando.regs[[x]])) {
    return(NULL)
  }
  split(pando.regs[[x]], f = pando.regs[[x]]$tf) %>% sapply(., "[[", "regions") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
dfRegionsPerReg.pando <- data.frame(
  No.regions.per.regulon = regionsPerReg.pando, 
  Method = rep("Pando", length(regionsPerReg.pando) )
)
rownames(dfGenesPerReg.pando) <- NULL
dfGenesPerReg.pando


# scMEGA
dfRegionsPerReg.scmega <- data.frame(
  No.regions.per.regulon = NA,
  Method = "scMEGA"
)


# Draw box plots to showcase No. of regions per regulon
nRegionsPerReg.df <- rbind(
  dfRegionsPerReg.stream, 
  dfRegionsPerReg.scplus, 
  dfRegionsPerReg.scenic, 
  dfRegionsPerReg.scglue, 
  dfRegionsPerReg.dinet, 
  dfRegionsPerReg.pando, 
  dfRegionsPerReg.scmega
)
nRegionsPerReg.df$Method <- factor(nRegionsPerReg.df$Method, 
                                 levels = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"))
p.nRegionsPerReg.df <- ggplot(nRegionsPerReg.df, aes(x = Method, y = No.regions.per.regulon, fill = Method)) + 
  geom_boxplot(lwd = 0.1) + scale_fill_manual(values = col.ls) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) + 
  theme(
    #axis.title = element_text(size = 24, color = "black"),        # Font settings for axis titles
    axis.text = element_text(size = 12, color = "black"),       # Font settings for axis text
    plot.title = element_text(size = 16, color = "black"), # Font settings for the plot title
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    panel.background = element_rect(fill = "white", linewidth = 0.3),
    legend.position = "none"
  )
nRegionsPerReg.median <- split(nRegionsPerReg.df, f = nRegionsPerReg.df$No.regions.per.regulon) %>% 
  sapply(., "[[", "No.regions.per.regulon") %>% sapply(., median) %>% round
qs::qsave(p.nRegionsPerReg.df, paste0(parent.dir, "boxplots_no_regions_per_regulon.qsave"))
png(filename = paste0(parent.dir, "boxplots_no_regions_per_regulon.png"), width = 300, height = 300, res = 150)
p.nRegionsPerReg.df
dev.off()


################################################################
#                                                              #
#                    Benchmark TF coverage                     #
#                                                              #
################################################################


benchTF.vis <- benchTF.sel[vis.data]
TfCov.ll <- list()
for (i in seq_along(vis.data)) {
  message (i)
  df <- data.frame(
    Method = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"),
    f.score = benchTF.vis[[i]][c(1:4, 7, 5, 6), 3]
  )
  df$Method <- factor(df$Method, levels = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"))
  p <- ggplot(df, aes(x = Method, y = f.score, fill = Method)) + 
    geom_bar(stat = "identity", fill = col.ls, lwd = 0.1) + scale_fill_manual(values = col.ls) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) + 
    theme(
      #axis.title = element_text(size = 24, color = "black"),        # Font settings for axis titles
      axis.text = element_text(size = 12, color = "black"),       # Font settings for axis text
      plot.title = element_text(size = 16, color = "black"), # Font settings for the plot title
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
      panel.background = element_rect(fill = "white", linewidth = 0.3),
      legend.position = "none"
    )
  TfCov.ll[[i]] <- p
}
library(ggpubr)
p.tf.cov <- ggarrange(plotlist = TfCov.ll, nrow = 2, ncol = 3)
qs::qsave(p.tf.cov, paste0(parent.dir, "barplots_benchmark_TF.qsave"))
png(filename = paste0(parent.dir, "barplots_benchmark_TF.png"), width = 1200, height = 900, res = 150)
p.tf.cov
dev.off()


################################################################
#                                                              #
#               Benchmark TF-to-region linkages                #
#                                                              #
################################################################


name.hash <- setNames(
  c("STREAM", "SCENIC+", "SCENIC", "GLUE", "Pando", "scMEGA", "DIRECT-NET"), 
  c("stream", "scenicplus", "scenic", "scglue", "pando", "scmega", "directnet")
)
benchTfRegion.vis <- benchTFRegion.final[vis.data]
TfRegion.ll <- list()
for (i in seq_along(vis.data)) {
  message (i)
  df <- data.frame(
    Method = name.hash[benchTfRegion.vis[[i]]$method],
    precision = as.numeric(benchTfRegion.vis[[i]][, 1])
  )
  df <- rbind(df, 
              c("SCENIC", 0), 
              c("DIRECT-NET", 0), 
              c("scMEGA", 0))
  df$Method <- factor(df$Method, levels = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"))
  df$precision <- as.numeric(df$precision)
  p <- ggplot(df, aes(x = Method, y = precision, fill = Method)) + 
    geom_boxplot(lwd = 0.1) + scale_fill_manual(values = col.ls) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) + 
    theme(
      #axis.title = element_text(size = 24, color = "black"),        # Font settings for axis titles
      axis.text = element_text(size = 12, color = "black"),       # Font settings for axis text
      plot.title = element_text(size = 16, color = "black"), # Font settings for the plot title
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
      panel.background = element_rect(fill = "white", linewidth = 0.3),
      legend.position = "none"
    )
  TfRegion.ll[[i]] <- p
}
p.tf.region <- ggarrange(plotlist = TfRegion.ll, nrow = 2, ncol = 3)
qs::qsave(p.tf.region, paste0(parent.dir, "barplots_benchmark_TF_region.qsave"))
png(filename = paste0(parent.dir, "barplots_benchmark_TF_region.png"), width = 1200, height = 900, res = 150)
p.tf.region
dev.off()


################################################################
#                                                              #
#               Benchmark region-to-gene linkages              #
#                                                              #
################################################################


benchRegionGene.vis <- benchRegionGene.sel[vis.data]
regionGene.ll <- list()
for (i in seq_along(vis.data)) {
  message (i)
  df <- data.frame(
    Method = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"),
    f.score = c(benchRegionGene.vis[[i]][1:2, 3], 0, benchRegionGene.vis[[i]][3, 3], 0, benchRegionGene.vis[[i]][4, 3], 0)
  )
  df$Method <- factor(df$Method, levels = c("STREAM", "SCENIC+", "SCENIC", "GLUE", "DIRECT-NET", "Pando", "scMEGA"))
  df$f.score <- as.numeric(df$f.score)
  p <- ggplot(df, aes(x = Method, y = f.score)) + 
    geom_bar(stat = "identity", fill = col.ls, lwd = 0.1) + scale_fill_manual(values = col.ls) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) + 
    theme(
      axis.text = element_text(size = 12, color = "black"),       # Font settings for axis text
      plot.title = element_text(size = 16, color = "black"), # Font settings for the plot title
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
      panel.background = element_rect(fill = "white", linewidth = 0.3),
      legend.position = "none"
    )
  regionGene.ll[[i]] <- p
}
p.region.gene <- ggarrange(plotlist = regionGene.ll, nrow = 2, ncol = 3)
qs::qsave(p.region.gene, paste0(parent.dir, "barplots_benchmark_region_gene.qsave"))
png(filename = paste0(parent.dir, "barplots_benchmark_region_gene.png"), width = 1200, height = 900, res = 150)
p.region.gene
dev.off()


################################################################
#                                                              #
#                      Panel arrangement                       #
#                                                              #
################################################################


png(filename = paste0(parent.dir, "benchmark_plots.png"), width = 2400, height = 1800, res = 150)
ggpubr::ggarrange(
  ggpubr::ggarrange(
    NA,
    ggpubr::ggarrange(p.nTfs, p.nGenesPerReg.df, p.nRegionsPerReg.df, ncol = 1, labels = LETTERS[2:4]),
    ggpubr::ggarrange(plotlist = p.tf.cov, labels = LETTERS[5:10]),
    nrow = 1, widths = c(1, 1, 2)
  ),
  ggpubr::ggarrange(
    ggpubr::ggarrange(plotlist = p.tf.region, labels = LETTERS[11:16]),
    ggpubr::ggarrange(plotlist = p.region.gene, labels = LETTERS[16:21]), 
    nrow = 1, widths = c(1, 1)
  ),
  ncol = 1
)
dev.off()


################################################################
#                                                              #
#                      Information of data                     #
#                                                              #
################################################################


save.image("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/5-benchmark-visualization.RData")
# load("/fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/5-benchmark-visualization.RData")