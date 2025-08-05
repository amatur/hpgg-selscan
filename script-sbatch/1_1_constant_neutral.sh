
#!/bin/bash
#SBATCH --job-name=1_1_constant_neutral       # Job name
#SBATCH --output=../1_constant/%x.out        # Output file (%j will be replaced with the job ID)
#SBATCH --error=../1_constant/%x.err          # Error file (%j will be replaced with the job ID)
#SBATCH --ntasks=1                    # Number of tasks (processes)
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=8G                      # Memory per node (10GB in this case)
#SBATCH --time=24:00:00               # Time limit (hh:mm:ss)
#SBATCH --mail-user=aur1111@psu.edu  # Email address for notifications
#SBATCH --mail-type=END,FAIL               # Send email on all events (BEGIN, END, FAIL, REQUEUE, etc.)

#SBATCH --account="zps5164_cr_default"          # Partition (queue) to submit to
#SBATCH --partition=basic

###SBATCH --account="zps5164_sc_default"          # Partition (queue) to submit to

source $(conda info --base)/etc/profile.d/conda.sh
conda activate test-slim

echo "Starting job at $(date)"  # Print the current date
cd ../1_constant
mkdir -p ../output/1_1_constant_neutral
slim 1_constant_neutral.slim
echo "Ending job at $(date)" 