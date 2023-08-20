#!/usr/bin/bash

DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/10X/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave
NAME=10X_hg38_Human_bone_marrow_with_T_cell_depleted
sbatch --job-name=STREAM_${NAME} --output=${NAME}.STREAM_out --export=DATA=$DATA 1-stream.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/scCATseq/scCATseq_hg19_K562_HCT116_HelaS3_PDX.qsave
NAME=scCATseq_hg19_K562_HCT116_HelaS3_PDX
sbatch --job-name=STREAM_${NAME} --output=${NAME}.STREAM_out --export=DATA=$DATA 1-stream.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/SHAREseq/SHAREseq_hg19_GM12878_rep2.qsave
NAME=SHAREseq_hg19_GM12878_rep2
sbatch --job-name=STREAM_${NAME} --output=${NAME}.STREAM_out --export=DATA=$DATA 1-stream.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/SHAREseq/SHAREseq_hg19_GM12878_rep3.qsave
NAME=SHAREseq_hg19_GM12878_rep3
sbatch --job-name=STREAM_${NAME} --output=${NAME}.STREAM_out --export=DATA=$DATA 1-stream.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/SNAREseq2/SNAREseq2_hg38_A549.qsave
NAME=SNAREseq2_hg38_A549
sbatch --job-name=STREAM_${NAME} --output=${NAME}.STREAM_out --export=DATA=$DATA 1-stream.sh
echo $NAME


DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/SNAREseq2/SNAREseq2_hg38_GM12878.qsave
NAME=SNAREseq2_hg38_GM12878
sbatch --job-name=STREAM_${NAME} --output=${NAME}.STREAM_out --export=DATA=$DATA 1-stream.sh
echo $NAME
