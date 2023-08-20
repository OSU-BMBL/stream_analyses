#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=46:58:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB

cd /fs/ess/PCON0022/guoqi/Yang/Stream/Directnet/output


DATA=/fs/scratch/PCON0022/liyang/GEX-ATAC-downsampled-benchmarking/10X/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave
NAME=10X_hg38_Human_bone_marrow_with_T_cell_depleted
sbatch --job-name=${NAME} --output=./log/$NAME.Direct_out --export=DATA=$DATA 1-run-directnet.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/GEX-ATAC-downsampled-benchmarking/SNAREseq2/SNAREseq2_hg38_A549.qsave
NAME=SNAREseq2_hg38_A549
sbatch --job-name=${NAME} --output=./log/$NAME.Direct_out --export=DATA=$DATA 1-run-directnet.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/GEX-ATAC-downsampled-benchmarking/SNAREseq2/SNAREseq2_hg38_GM12878.qsave
NAME=SNAREseq2_hg38_GM12878
sbatch --job-name=${NAME} --output=./log/$NAME.Direct_out --export=DATA=$DATA 1-run-directnet.sh
echo $NAME
