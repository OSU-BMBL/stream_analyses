#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60GB


source ~/.bashrc
module load python/3.9-2022.05
cd /fs/scratch/PCON0022/liyang/STREAM-revision/Benchmarking/scglue
code=/fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/scglue


conda activate /fs/ess/PCON0022/liyang/.conda/scglue


python ${code}/1-run-scglue.py $DATA
sleep 0.1s
echo "End running scglue on ${DATA}"
echo ""
echo ""


conda deactivate
