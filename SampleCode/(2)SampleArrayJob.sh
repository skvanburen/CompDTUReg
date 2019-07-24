#!/bin/bash



#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g

mkdir -p SaveInfRepsAsRDataLogs
srun R CMD BATCH (2)SaveInfRepsAsRData.R SaveInfRepsAsRDataLogs/RLog_$SLURM_ARRAY_TASK_ID.Rout
