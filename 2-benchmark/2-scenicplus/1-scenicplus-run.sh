#!/usr/bin/bash


cd /fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/scenicplus
while read DATA
do
    echo $DATA
    sbatch --job-name=scenicplus_${DATA} --output=${DATA}.scenicplus_out --export=DATA=$DATA 1-scenicplus.sh
    sleep 0.1s
done < dataset_list.txt
