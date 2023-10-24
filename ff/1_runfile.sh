#!/bin/bash
#SBATCH --job-name=all
#SBATCH --nodes=1
#SBATCH --time=0-36
#SBATCH --account=isda
#SBATCH --mail-user=emmanuel.branlard@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails
#SBATCH -o slurm-%x.log                      # Output %j: job number
echo "# Working directory:" $SLURM_SUBMIT_DIR
echo "# Job name:" $SLURM_JOB_NAME
echo "# Job ID: " $SLURM_JOBID
echo "# Starting job at: " $(date)

export OMP_NUM_THREADS=26

export OMP_PROC_BIND=spread
export KMP_AFFINITY=balanced

#ffbin='/home/ebranlar/_bin/FAST.Farm-wat-gcc'
#source /home/ebranlar/_env/ebra-gcc.sh

ffbin='/home/ebranlar/_bin/FAST.Farm-wat-intel'
source /home/ebranlar/_env/ebra-intel.sh

echo "# Binary: " $ffbin

cd neutral
$ffbin FF-WAT.fstf 2>&1  &
sleep 5
$ffbin FF-NoWAT.fstf 2>&1  &
cd ..

cd stable
$ffbin FF-WAT.fstf 2>&1  &
sleep 5
$ffbin FF-NoWAT.fstf 2>&1  &
cd ..

wait


echo "# Ending job at: " $(date)
