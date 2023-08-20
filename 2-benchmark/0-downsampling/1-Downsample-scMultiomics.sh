#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=10:29:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB
#SBATCH --job-name=1-Downsample-scMultiomics
#SBATCH --output=1-Downsample-scMultiomics.txt


module load R/4.2.1-gnu11.2
cd /fs/ess/PCON0022/liyang/STREAM-revision/Downsampling


# Start the timer
start_time=$(date +%s)


Rscript 1-Downsample-scMultiomics.R


# End the timer
end_time=$(date +%s)


# Calculate the running time
running_time=$((end_time - start_time))


echo "Script execution time: $running_time seconds"