#####################################################
#                                                   #
# Calculate pairwise similarity between eRegulons   #
# in different time points                          #
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



# Load files
obj <- qs::qread(paste0(R.dir, "Object_Shane_annot.qsave"))
dim(obj)
table(obj$celltypes)
en.reg.2.5 <- qs::qread(paste0(R.dir, "Extended_eGRNs_2.5.qsave"))
en.reg.5.7 <- qs::qread(paste0(R.dir, "Extended_eGRNs_5.7.qsave"))
en.reg.13 <- qs::qread(paste0(R.dir, "Extended_eGRNs_13.qsave"))
length(en.reg.2.5)
length(en.reg.5.7)
length(en.reg.13)



# Build eRegulon similarity between 2.5 and 5.7 months
sim.df1 <- data.frame(matrix(ncol = 5, nrow = 0)) # similarity data frame between 2.5 and 5.7 months
for (i in seq_along(en.reg.2.5)) {
  i.links <- paste(GRangesToString(en.reg.2.5[[i]]$links), 
                   en.reg.2.5[[i]]$links$gene, sep = "_")
  j.genes <- unique(en.reg.2.5[[i]]$links$gene)
  i.enhs <- unique(GRangesToString(en.reg.2.5[[i]]$links))
  for (j in seq_along(en.reg.5.7)) {
    j.links <- paste(GRangesToString(en.reg.5.7[[j]]$links), 
                     en.reg.5.7[[j]]$links$gene, sep = "_")
    j.genes <- unique(en.reg.5.7[[j]]$links$gene)
    j.enhs <- unique(GRangesToString(en.reg.5.7[[j]]$links))
    
    ratio.links <- length(intersect(i.links, j.links)) / min(length(i.links), length(j.links))
    ratio.genes <- length(intersect(i.genes, j.genes)) / min(length(i.genes), length(j.genes))
    ratio.enhs <- length(intersect(i.enhs, j.enhs)) / min(length(i.enhs), length(j.enhs))
    
    sim.df1 <- rbind(sim.df1, c(i, j, ratio.links, ratio.genes, ratio.enhs))
  }
}
colnames(sim.df1) <- c("id1", "id2", "link", "gene", "enhancer")
sim.df1$id1 <- paste0("eR_2.5_", sim.df1$id1)
sim.df1$id2 <- paste0("eR_5.7_", sim.df1$id2)
dim(sim.df1)
head(sim.df1)
apply(sim.df1, 2, range)
qs::qsave(sim.df1, paste0(R.dir, "similarity_2.5_5.7_mon.qsave"))



# Build eRegulon similarity between 5.7 months and 13+ months
sim.df2 <- data.frame(matrix(ncol = 5, nrow = 0)) # similarity data frame between 5.7 and 13+ months
for (i in seq_along(en.reg.5.7)) {
  i.links <- paste(GRangesToString(en.reg.5.7[[i]]$links), 
                   en.reg.5.7[[i]]$links$gene, sep = "_")
  j.genes <- unique(en.reg.5.7[[i]]$links$gene)
  i.enhs <- unique(GRangesToString(en.reg.5.7[[i]]$links))
  for (j in seq_along(en.reg.13)) {
    j.links <- paste(GRangesToString(en.reg.13[[j]]$links), 
                     en.reg.13[[j]]$links$gene, sep = "_")
    j.genes <- unique(en.reg.13[[j]]$links$gene)
    j.enhs <- unique(GRangesToString(en.reg.13[[j]]$links))
    
    ratio.links <- length(intersect(i.links, j.links)) / min(length(i.links), length(j.links))
    ratio.genes <- length(intersect(i.genes, j.genes)) / min(length(i.genes), length(j.genes))
    ratio.enhs <- length(intersect(i.enhs, j.enhs)) / min(length(i.enhs), length(j.enhs))
    
    sim.df2 <- rbind(sim.df2, c(i, j, ratio.links, ratio.genes, ratio.enhs))
  }
}
colnames(sim.df2) <- c("id1", "id2", "link", "gene", "enhancer")
sim.df2$id1 <- paste0("eR_5.7_", sim.df2$id1)
sim.df2$id2 <- paste0("eR_13_", sim.df2$id2)
dim(sim.df2)
head(sim.df2)
apply(sim.df2, 2, range)
qs::qsave(sim.df2, paste0(R.dir, "similarity_5.7_13_mon.qsave"))



# Convert data frame to igraph and Cytoscape object
sim.df <- rbind(sim.df1, sim.df2)
ig.lst <- lapply(3 : 5, function(i) {
  tmp.df <- sim.df[, c(1:2, i)]
  tmp.ig <- graph_from_data_frame(sim.df, directed = F)
  E(tmp.ig)$weight <- sim.df[, i]
  vcount(tmp.ig)
  ecount(tmp.ig)
  tmp.ig
})
names(ig.lst) <- c("link", "gene", "enhancer")
E(ig.lst[[1]])$weight %>% head
head(sim.df)
qs::qsave(ig.lst, paste0(R.dir, "igraph_list_link_gene_enh.qsave"))



# Plot igraph with linkage similarity as weights (n = 50))
tf.lst <- c(
  "Runx1",  
  "Fos", 
  "Nfe2", 
  "Jund", 
  "Fosl2", 
  "Nr2c2", 
  "Esr2", 
  "Nkx3-1", 
  "Elk4", 
  "Ar", 
  "Tlx1", 
  "Jun", 
  "Rora", 
  "Zbtb33", 
  "Plag1"
)

sim.link <- sim.df[order(sim.df$link, decreasing = T), 1:3]
sim.link
sim.link$id1 <- paste0("d", strsplit(sim.link$id1, split = "_") %>% sapply(., "[[", 2), "_", 
                       strsplit(sim.link$id1, split = "_") %>% sapply(., function(x) {
  x[3] <- get(paste0("en.reg.", x[2]))[as.numeric(x[3])][[1]]$TF %>% 
    strsplit(., split = "::") %>% sapply(., function(x) {tail(x, n = 1)}) %>% 
    tolower %>% capitalize
}))
sim.link$id2 <- paste0("d", strsplit(sim.link$id2, split = "_") %>% sapply(., "[[", 2), "_", 
  strsplit(sim.link$id2, split = "_") %>% sapply(., function(x) {
  x[3] <- get(paste0("en.reg.", x[2]))[as.numeric(x[3])][[1]]$TF %>% 
    strsplit(., split = "::") %>% sapply(., function(x) {tail(x, n = 1)}) %>% 
    tolower %>% capitalize
}))
sim.link$link <- as.numeric(sim.link$link)
rownames(sim.link) <- NULL
dim(sim.link)
sim.link <- sim.link[apply(sim.link, 1, function(x) {
  ar <- strsplit(unlist(x[1:2]), split = "_")
  if (ar[[1]][2] %in% tf.lst & ar[[2]][2] %in% tf.lst & (x[3] > 0)) {
    return(T)
  } else {
    return(F)
  }
}),]
dim(sim.link)
sim.link <- sim.link[order(sim.link$link, decreasing = T),]
sim.link <- sim.link[sim.link$link > 0,]
dim(sim.link)
head(sim.link)
tail(sim.link)
top.sim <- sim.link[1:556,]
top.sim
ig.df <- distinct(top.sim, id1, id2, .keep_all = T)
dim(ig.df)
setdiff(tf.lst, strsplit(c(ig.df$id1, ig.df$id2), split = "_") %>% 
          sapply(., "[[", 2) %>% unique)
tf.lst
ig.link <- graph_from_data_frame(ig.df, directed = T)
E(ig.link)$weight <- ig.df$link
head(E(ig.link)$weight)
head(ig.df)
tail(ig.df)
range(ig.df$link)
vcount(ig.link)
ecount(ig.link)
layer <- rep(1, vcount(ig.link))
layer[grep("^d5.7_", V(ig.link)$name)] <- 2
layer[grep("^d13_", V(ig.link)$name)] <- 3
layout <- layout_with_sugiyama(ig.link, layers = layer)
plot(ig.link,
     layout=cbind(layer, layout$layout[,1]),
     vertex.shape=c("circle","circle","circle")[layer],
     vertex.size=c(50,20,0)[layer],
     vertex.label.dist=c(0,0,.8)[layer],
     vertex.label.degree=0)
qs::qsave(ig.link, paste0(R.dir, "en_reg_sim.qsave"))


save.image(paste0(scratch.dir, "21_eReg_similarity.RData"))
