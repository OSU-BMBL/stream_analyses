######################################################
#                                                    #
#   Extract the gene symbols and save them in txt    #
#                                                    #
######################################################


# Parameters
data_dir = "/fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/10X-Genomics-Multiome/10X-human-bone-marrow/"


# Modules
import scanpy as sc


# Load annData
adata = sc.read_h5ad(data_dir + "bm_multiome_rna.h5ad")


# Gene symbol conversion
gene_names = adata.var_names.tolist()

with open(data_dir + 'bm_multiome_rna.txt', 'w') as f:
    for item in gene_names:
        f.write("%s\n" % item)
