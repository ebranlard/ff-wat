#!/bin/bash
#SBATCH --job-name=neut-ff
#SBATCH --nodes=1
#SBATCH --time=0-01
#SBATCH --account=isda
#SBATCH --mail-user=emmanuel.branlard@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails
#SBATCH -o slurm-%x-%j.log                      # Output

cores=$(( 36*$SLURM_JOB_NUM_NODES )); fi

echo "# Working directory:" $SLURM_SUBMIT_DIR
echo "# Job name:" $SLURM_JOB_NAME
echo "# Job ID: " $SLURM_JOBID
echo "# Submit time is" $(squeue -u $USER -o '%30j %20V' | grep -e $SLURM_JOB_NAME | awk '{print $2}')
echo "# Starting job at: " $(date)
echo "# Using: " $cores "core(s)"

module purge
module purge
module load comp-intel mkl


ffbin= '/home/ebranlar/_bin/FAST.Farm-v3.4-intel-omp'
input= FF-NoWAT.fstf

$ffbin $input 2>&1

echo "# Ending job at: " $(date)
