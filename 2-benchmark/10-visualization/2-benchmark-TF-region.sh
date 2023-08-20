#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=06:29:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50GB
#SBATCH --job-name=2-benchmark-TF-region
#SBATCH --output=2-benchmark-TF-region.txt


cd /fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/
module load R/4.2.1-gnu11.2


# Start the timer
start_time=$(date +%s)
Rscript	2-benchmark-TF-region.R


# End the timer
end_time=$(date +%s)


# Calculate the running time
running_time=$((end_time - start_time))


# Display the running time
echo "Script execution time: $running_time seconds"
