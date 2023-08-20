#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=46:58:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB


module load R/4.2.1-gnu11.2
cd /fs/ess/PCON0022/liyang/STREAM-revision/Benchmarking/stream/


# Start the timer
start_time=$(date +%s)


# Run STREAM
Rscript 1-stream.R $DATA


# End the timer
end_time=$(date +%s)


# Calculate the running time
running_time=$((end_time - start_time))


# Display the running time
echo "Script execution time: $running_time seconds"