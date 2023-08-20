######################################################
#                                                    #
#             Analysis on multiple genes             #
#                                                    #
######################################################



# Rdata
load("6_DSLL_single_gene_analyses.RData")



# Candidate TF LHX2
normal.tfs <- normal.en.grns$TF %>% unique
prol.tfs <- prol.en.grns$TF %>% unique
tumor.tfs <- tumor.en.grns$TF %>% unique
Reduce("intersect", list(normal.tfs, prol.tfs, tumor.tfs))



# Get diff genes
sig.dars <- de.genes[de.genes$avg_log2FC > 0 & de.genes$p_val_adj < 0.05,]
prol.genes <- intersect(prol.en.grns[prol.en.grns$TF == "LHX2", "gene"],
                        unique(sig.dars[sig.dars$cluster == "Tumor B prol", "gene"]))
tumor.genes <- intersect(tumor.en.grns[tumor.en.grns$TF == "LHX2", "gene"],
                        unique(sig.dars[sig.dars$cluster == "Tumor B", "gene"]))
normal.genes <- setdiff(normal.en.grns[normal.en.grns$TF == "LHX2", "gene"],
                        union(prol.genes, tumor.genes))



# Get the candidate enhs
normal.enhs <- normal.en.grns[normal.en.grns$TF == "LHX2", "enhancer"]
prol.enhs <- prol.en.grns[prol.en.grns$TF == "LHX2", "enhancer"]
tumor.enhs <- tumor.en.grns[tumor.en.grns$TF == "LHX2", "enhancer"]
join.enhs <- Reduce("intersect", list(normal.enhs, prol.enhs, tumor.enhs))



join.da.enhs <- b.da.enhs[b.da.enhs$gene %in% join.enhs & 
                            b.da.enhs$avg_log2FC > 0 & 
                            b.da.enhs$p_val_adj < 0.05,]
key.enh <- join.da.enhs$gene[1]
da.enhs[da.enhs$gene %in% join.enhs & 
          da.enhs$avg_log2FC > 0 & 
          da.enhs$p_val_adj < 0.05,]
normal.diffGenes <- normal.en.grns[normal.en.grns$enhancer == key.enh & 
                                 normal.en.grns$TF == "LHX2", "gene"] %>% unique
prol.diffGenes <- prol.en.grns[prol.en.grns$enhancer == key.enh & 
                             prol.en.grns$TF == "LHX2", "gene"] %>% unique
tumor.diffGenes <- tumor.en.grns[tumor.en.grns$enhancer == key.enh & 
                               tumor.en.grns$TF == "LHX2", "gene"] %>% unique
intersect(prol.diffGenes, unique(b.de.genes$gene))


# Prepare TF binding sites on the key enhancer
load("/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/TFBS_list.rda")
key.tf <- "LHX2"
require(Signac)
key.tf.sites <- pblapply(list(key.tf), function(x) {
  tf.sites <- TFBS.list$Human$peak[which(TFBS.list$Human$TF == x)]
  subject.query <- findOverlaps(StringToGRanges(key.enh), tf.sites)
  key.sites <- tf.sites[subjectHits(subject.query) %>% unique]
})[[1]]
DefaultAssay(b.obj) <- "ATAC"
b.colors <- c("#FFADAD", "#FF5D5D", "#D61C4E")
show_col(b.colors)
names(b.colors) <- unique(b.obj$ct.label) %>% rev
show_col(b.colors)
tf.colors <- "#3F52E3"
mcols(key.tf.sites)$color <- rep(tf.colors, length(key.tf.sites))



# Get enh-gene linkages
gene.regions <- obj@assays$ATAC@annotation
key.gene <- tumor.diffGenes
gene.regions <- get_gene_loc(key.gene = key.gene, gene.regions = gene.regions)
links <- punion(resize(gene.regions, 1, "start"), rep(resize(StringToGRanges(key.enh), 1, "center"), 
                                  length(gene.regions)), 
                     fill.gap = T, ignore.strand = T)
key.region <- paste0("chr3-", min(start(links)), "-", 
                     max(end(links))) %>% StringToGRanges
region.highlight <- resize(key.tf.sites, width = 3000)
mcols(links)$gene <- gene.regions$gene_name
rna.m <- GetAssayData(b.obj, assay = "RNA", slot = "data")
atac.m <- GetAssayData(b.obj, assay = "ATAC", slot = "data")
link.scores <- rbindlist(lapply(rev(unique(b.obj$ct.label)), function(x) {
  enh.cells <- names(which((atac.m > 0)[key.enh, b.obj@meta.data[, "ct.label"] == x]))
  genes <- lapply(key.gene, function(y) {
    yy <- names(which((rna.m > 0)[y, b.obj@meta.data[, "ct.label"] == x]))
    length(intersect(enh.cells, yy)) / length(enh.cells)
  })
  names(genes) <- key.gene
  genes
}))
rownames(link.scores) <- levels(b.obj$ct.label)[7:9]


mcols(links)$score <- link.scores[1,] %>% unlist
mcols(links)$peak <- rep(GRangesToString(region.highlight), length(gene.regions))
DefaultAssay(b.obj) <- "ATAC"
Links(b.obj) <- links
p.multiGeneCov.normal <- get_coverage_plot(object = b.obj, 
                                     region = key.region, 
                                     features = key.gene, links = T,
                                     ranges.group.by = "ct.label", peaks = T, 
                                     bigwig = NULL, region.highlight = region.highlight, 
                                     colors = b.colors, size = 10)
qs::qsave(p.multiGeneCov.normal, paste0(R.dir, "Coverage_plot_multi_genenormal.qsave"))
save_image(p = p.multiGeneCov.normal, path = paste0(image.dir, "Coverage_plot_multi_genenormal.png"), 
           width = 3600, height = 1500)


mcols(links)$score <- link.scores[2,] %>% unlist
mcols(links)$peak <- rep(GRangesToString(region.highlight), length(gene.regions))
prol.obj <- b.obj
DefaultAssay(prol.obj) <- "ATAC"
Links(prol.obj) <- links
p.multiGeneCov.prol <- get_coverage_plot(object = prol.obj, 
                                    region = key.region, 
                                    features = key.gene, links = T,
                                    ranges.group.by = "ct.label", peaks = T, 
                                    bigwig = NULL, region.highlight = region.highlight, 
                                    colors = b.colors, size = 10)
qs::qsave(p.multiGeneCov.prol, paste0(R.dir, "Coverage_plot_multi_gene_prol.qsave"))
save_image(p = p.multiGeneCov.prol, path = paste0(image.dir, "Coverage_plot_multi_gene_prol.png"), 
           width = 3600, height = 1500)


mcols(links)$score <- link.scores[3,] %>% unlist
mcols(links)$peak <- rep(GRangesToString(region.highlight), length(gene.regions))
tumor.obj <- b.obj
DefaultAssay(tumor.obj) <- "ATAC"
Links(tumor.obj) <- links
p.multiGeneCov.tumor <- get_coverage_plot(object = tumor.obj, 
                                         region = key.region, 
                                         features = key.gene, links = T,
                                         ranges.group.by = "ct.label", peaks = T, 
                                         bigwig = NULL, region.highlight = region.highlight, 
                                         colors = b.colors, size = 10)
qs::qsave(p.multiGeneCov.tumor, paste0(R.dir, "Coverage_plot_multi_gene_tumor.qsave"))
save_image(p = p.multiGeneCov.tumor, path = paste0(image.dir, "Coverage_plot_multi_gene_tumor.png"), 
           width = 3600, height = 1500)



save.image("/fs/ess/scratch/PCON0022/liyang/stream/Case_2_AD/7_DSLL_multi_gene_analyses.Rdata")