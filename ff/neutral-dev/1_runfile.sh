#!/bin/bash
#SBATCH --job-name=neut-db
#SBATCH --nodes=1
#SBATCH --time=0-00:10
#SBATCH --account=isda
#SBATCH --mail-user=emmanuel.branlard@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails
#SBATCH -o slurm-%x.log                      # Output
echo "# Working directory:" $SLURM_SUBMIT_DIR
echo "# Job name:" $SLURM_JOB_NAME
echo "# Job ID: " $SLURM_JOBID
echo "# Starting job at: " $(date)

export OMP_NUM_THREADS=10

#ffbin='/home/ebranlar/_bin/FAST.Farm-dev-gcc'
#source /home/ebranlar/_env/ebra-gcc.sh

ffbin='/home/ebranlar/_bin/FAST.Farm-dev-intel'
source /home/ebranlar/_env/ebra-intel.sh

echo "# Binary: " $(ffbin)

$ffbin FF.fstf 2>&1

echo "# Ending job at: " $(date)
