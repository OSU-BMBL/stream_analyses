#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=46:58:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB

echo $NAME
module load R/4.1.0-gnu9.1
cd /fs/ess/PCON0022/guoqi/Yang/Stream/Directnet/output
#DATA=/fs/scratch/PCON0022/liyang/STREAM-revision/Downsampled-data/10X/10X_mm10_Mouse_kidney_nuclei_isolated.qsave 

# Start the timer
start_time=$(date +%s)


# Run Directnet
Rscript directnet.r $DATA


# End the timer
end_time=$(date +%s)


# Calculate the running time
running_time=$((end_time - start_time))


# Display the running time
echo "Script execution time: $running_time seconds"