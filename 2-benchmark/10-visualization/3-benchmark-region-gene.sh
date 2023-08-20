#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=11:29:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50GB
#SBATCH --job-name=3-benchmark-region-gene
#SBATCH --output=3-benchmark-region-gene.txt


cd /fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/
module load R/4.2.1-gnu11.2


# Start the timer
start_time=$(date +%s)
Rscript	3-benchmark-region-gene.R


# End the timer
end_time=$(date +%s)


# Calculate the running time
running_time=$((end_time - start_time))


# Display the running time
echo "Script execution time: $running_time seconds"
