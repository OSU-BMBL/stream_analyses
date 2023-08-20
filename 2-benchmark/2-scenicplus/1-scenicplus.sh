#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50GB


source ~/.bashrc
module load python/3.9-2022.05
cd /fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scenicplus
code=/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/scenicplus


conda activate scenicplus


python ${code}/1-run-scenicplus-minimal.py $DATA


conda deactivate
