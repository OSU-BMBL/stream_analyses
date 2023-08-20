#!/usr/bin/bash


cd /fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/scglue
while read DATA
do
    echo $DATA
    sbatch --job-name=scglue_${DATA} --output=${DATA}.scglue_out --export=DATA=$DATA 1-scglue.sh
    sleep 0.1s
done < dataset_list.txt
