#!/bin/bash
#SBATCH --job-name=stab-wa
#SBATCH --nodes=1
#SBATCH --time=0-36
#SBATCH --account=isda
#SBATCH --mail-user=emmanuel.branlard@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails
#SBATCH -o slurm-%x.log                      # Output %j: job number

cores=$(( 36*$SLURM_JOB_NUM_NODES ))

echo "# Working directory:" $SLURM_SUBMIT_DIR
echo "# Job name:" $SLURM_JOB_NAME
echo "# Job ID: " $SLURM_JOBID
echo "# Starting job at: " $(date)
echo "# Using: " $cores "core(s)"

# ffbin='/home/ebranlar/_bin/FAST.Farm-vWAT'
# module purge
# module load craype-x86-spr
# module load intel-oneapi-mkl/2023.2.0-intel
# module load intel-oneapi-mpi/2021.10.0-intel
# module load intel-oneapi-compilers/2023.2.0
# module load hdf5/1.14.1-2-intel-oneapi-mpi-intel

ffbin='/home/ebranlar/_bin/FAST.Farm-vWAT-Regis-Nrel'
source /home/ebranlar/_env/regis-nrel.sh

$ffbin FF-WAT.fstf  2>&1


echo "# Ending job at: " $(date)