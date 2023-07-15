#!/bin/bash
#SBATCH --job-name=s8at150_0.1z0_2.5m
#SBATCH --output amr1run.log.%j
#SBATCH --nodes=43
#SBATCH --time=3-00
#SBATCH --account=car
#SBATCH --mail-user=registhedin@gmail.com
#SBATCH --mail-type=ALL

# This stable problem has 6144 grids
# 6144/(36*32) = 5.6 grids
# let's do 4 grids/code. 6144/4 = 1536 cores (43 nodes == 43*36=1548 cores)
cores=1536

source $HOME/.bash_profile
#cores=$SLURM_NTASKS
#if [ -z $cores ]; then cores=$(( 36*$SLURM_JOB_NUM_NODES )); fi

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

amrbin='/home/rthedin/repos/amr-wind/build/amr_wind'

export EXAWIND_DIR=/nopt/nrel/ecom/exawind/exawind-2020-09-21/install/gcc
export MPI_TYPE_DEPTH=15

#rm -rf post_processing 
srun -n $cores --cpu_bind=cores $amrbin setup_precursor_stable.i > log.amr_wind.abl 2>&1

echo "Ending AMR-wind job at: " $(date)