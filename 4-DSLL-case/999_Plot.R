# Arrange panels




# Parameters
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Images/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Rdata/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
source(paste0(tool.dir, "visual_tools.R"))




p.act.umap <- qs::qread(paste0(R.dir, "UMAP_DeepMAPS_active_celltypes.qsave"))
p.bar.ct.eRegs <- qs::qread(paste0(R.dir, "Barplot_celltype_eRegulons.qsave"))
p.en.reg.info <- qs::qread(paste0(R.dir, "Boxplots_eReghulon_basic_info.qsave"))
p.hist.gene.enh <- qs::qread(paste0(R.dir, "Histograms_genes_enhs.qsave"))
# Venn diagram
# LHX2 network x 2
p.vln.prol.lhx2 <- qs::qread(paste0(R.dir, "VlnPlot_prol_LHX2.qsave"))
# EGR1 network x 2
p.heatmap.prol.chipseq <- qs::qread(paste0(R.dir, "Heatmap_overlap_ChIP_seq.qsave")) %>% as.ggplot
p.heatmap.tumor.chipseq <- qs::qread(paste0(R.dir, "Heatmap_overlap_ChIP_seq.qsave")) %>% as.ggplot
require(ggpubr)
p <- ggarrange(
  ggarrange(
    ggarrange(p.act.umap, p.bar.ct.eRegs, 
            widths = c(2, 1), labels = LETTERS[1:2], nrow = 1), 
    ggarrange(p.en.reg.info, p.hist.gene.enh, ncol = 1, heights = c(1, 2), 
              labels = LETTERS[3:4]), nrow = 1, widths = c(3, 1)
    ),
  ggarrange(
    ggarrange(
    ggarrange(NA, NA, NA, 
              nrow = 1, labels = LETTERS[5:7]),
    p.vln.prol.lhx2,
    # ggarrange(p.vln.prol.lhx2, p.heatmap.prol.chipseq, widths = c(2, 1), nrow = 1), 
    labels = LETTERS[8], ncol = 1, heights = c(1, 2)
  ), 
  ggarrange(NA, NA, 
            ggarrange(p.heatmap.prol.chipseq, p.heatmap.tumor.chipseq, 
                      nrow = 1, labels = LETTERS[11:12]), ncol = 1, 
            labels = LETTERS[c(9:10, NA, NA)]), 
  widths = c(3, 2)
  ), 
  ncol = 1, heights = c(1, 2)
)


save_image(p = p, path = paste0(image.dir, "DSLL.png"), 
           width = 6000, height = 6000)
