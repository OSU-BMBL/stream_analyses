# Modules
import os, sys



# Source scripts
genome_funs = "/fs/ess/PCON0022/liyang/Python_utilities/Functions/"
sys.path.append(genome_funs)
import Genome_tools



# Parameters
bw_files = [
  "ENCFF042MFJ.bigWig", # SPI1
  "ENCFF011FHF.bigWig", # STAT3
  "ENCFF579VGX.bigWig" # TCF12
]


# Preparation
os.chdir("/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Tables/")
region_file = "/fs/ess/PCON0022/liyang/STREAM/Case_1_DSLL_cell_type/Tables/All_enhs.txt"


# Run Python 1
Genome_tools.calc_profile(region_file, bw_files, width = 750, val_type = "max", nBins = 50)






# Run Python 2
# The folowing codes aim to generate the TF binding profiles centered at the midpoint
# with diameter 750 bp
lymph_bw = [
  "lymph_node_lymphoma_14k_atac_cut_sites.bigWig"
]
Genome_tools.calc_profile(region_file, lymph_bw, width = 750, 
                          val_type = "max", nBins = 50)
