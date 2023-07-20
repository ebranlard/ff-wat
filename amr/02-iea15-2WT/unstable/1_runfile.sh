#!/bin/bash
#SBATCH --job-name=unstable_2wt
#SBATCH --nodes=72
#SBATCH --time=8-00
#SBATCH --account=mmc
#SBATCH --mail-user=emmanuel.branlard@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails
#SBATCH -o slurm-%x-%j.log                      # Output

# This problem has 12800 grids
# if 4 grids per core, we need 12800/5 = 2560 cores
# 72 nodes contain 2592 cores
#==> Ask for 72 nodes, use 2560 cores (72*36=2592)

source $HOME/.bash_profile
#cores=$SLURM_NTASKS
cores=2560
if [ -z $cores ]; then cores=$(( 36*$SLURM_JOB_NUM_NODES )); fi

echo "Working directory is" $SLURM_SUBMIT_DIR
echo "Job name is" $SLURM_JOB_NAME
echo "Job ID is " $SLURM_JOBID
echo "Submit time is" $(squeue -u $USER -o '%30j %20V' | grep -e $SLURM_JOB_NAME | awk '{print $2}')
echo "Starting AMR-wind job at: " $(date)
echo "using" $cores "core(s)"

module purge
module load gcc
module load mpt
module load cmake
module load mkl
module load netcdf-c/4.7.3

amrbin='/home/rthedin/repos/amr-wind_2023_05_01_openfastrestart_bugfix/amr-wind/build_main_2023_05_01_10cd537/amr_wind'
input=amr.i

export EXAWIND_DIR=/nopt/nrel/ecom/exawind/exawind-2020-09-21/install/gcc
export MPI_TYPE_DEPTH=15

# Ran this manually:
# python calc_inflowoutflow_stats.py -sf ../../02_precursor_shell/unstable.W.8at150.20dTinv_0.05q_0.75z0_850zi_3.84x1.28x0.9km_res2.5m/post_processing/abl_statistics76826.nc -ts 15000 -te 16800 -if setup_turbine_unstable.startAt15000.i

#rm -rf post_processing 
srun -n $cores --cpu_bind=cores $amrbin $input 2>&1

echo "Ending AMR-wind job at: " $(date)
