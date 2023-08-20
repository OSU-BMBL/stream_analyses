#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=11:29:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB
#SBATCH --job-name=1-Download-10X-human-bone-marrow
#SBATCH --output=1-Download-10X-human-bone-marrow.txt


module load R/4.2.1-gnu11.2
cd /fs/scratch/PCON0022/liyang/Joint-ATAC-RNA/10X-human-bone-marrow/


# Start the timer
start_time=$(date +%s)


# Rscript 1-Download-scNOMeRe-seq.R
wget https://zenodo.org/record/6383269/files/BM_CD34_Rep1_atac_fragments.tsv.gz
wget https://zenodo.org/record/6383269/files/BM_CD34_Rep1_atac_fragments.tsv.gz.tbi
wget https://zenodo.org/record/6383269/files/BM_CD34_Rep1_filtered_feature_bc_matrix.h5
wget https://zenodo.org/record/6383269/files/BM_CD34_Rep2_atac_fragments.tsv.gz
wget https://zenodo.org/record/6383269/files/BM_CD34_Rep2_atac_fragments.tsv.gz.tbi
wget https://zenodo.org/record/6383269/files/BM_CD34_Rep2_filtered_feature_bc_matrix.h5
wget https://zenodo.org/record/6383269/files/bm_multiome_atac.h5ad
wget https://zenodo.org/record/6383269/files/bm_multiome_rna.h5ad
wget https://zenodo.org/record/6383269/files/BM_Tcelldep_Rep1_atac_fragments.tsv.gz
wget https://zenodo.org/record/6383269/files/BM_Tcelldep_Rep1_atac_fragments.tsv.gz.tbi
wget https://zenodo.org/record/6383269/files/BM_Tcelldep_Rep1_filtered_feature_bc_matrix.h5
wget https://zenodo.org/record/6383269/files/BM_Tcelldep_Rep2_atac_fragments.tsv.gz
wget https://zenodo.org/record/6383269/files/BM_Tcelldep_Rep2_atac_fragments.tsv.gz.tbi
wget https://zenodo.org/record/6383269/files/BM_Tcelldep_Rep2_filtered_feature_bc_matrix.h5
wget https://zenodo.org/record/6383269/files/cd34_multiome_atac.h5ad
wget https://zenodo.org/record/6383269/files/cd34_multiome_rna.h5ad


# End the timer
end_time=$(date +%s)


# Calculate the running time
running_time=$((end_time - start_time))


# Display the running time
echo "Script execution time: $running_time seconds"