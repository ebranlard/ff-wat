#!/bin/bash
#SBATCH --job-name=neut_0wt
#SBATCH --nodes=56
#SBATCH --time=8-00
#SBATCH --account=isda
#SBATCH --mail-user=emmanuel.branlard@nrel.gov
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails
#SBATCH -o slurm-%x-%j.log                      # Output

# This problem has 7936 grids
# if 4 grids per core, we need 7936/4 = 1984 cores
# 56 nodes contain 2016 cores
#==> Ask for 56 nodes, use 1984 cores (56*36=2016)

#source $HOME/.bash_profile
#cores=$SLURM_NTASKS
cores=1984
if [ -z $cores ]; then cores=$(( 36*$SLURM_JOB_NUM_NODES )); fi

echo "Working directory is" $SLURM_SUBMIT_DIR
echo "Job name is" $SLURM_JOB_NAME
echo "Job ID is " $SLURM_JOBID
echo "Submit time is" $(squeue -u $USER -o '%30j %20V' | grep -e $SLURM_JOB_NAME | awk '{print $2}')
echo "Starting AMR-wind job at: " $(date)
echo "using" $cores "core(s)"

module purge
module load gcc/8.4.0
module load mpt
module load cmake
module load mkl
module load netcdf-c/4.7.3


#amrbin='/home/rthedin/repos/amr-wind/build_main_2022_02_02_ef466d9/amr_wind'
amrbin='/home/ebranlar/amr-wind-bin/amr_wind'
input=amr.i

export EXAWIND_DIR=/nopt/nrel/ecom/exawind/exawind-2020-09-21/install/gcc
export MPI_TYPE_DEPTH=15

# Ran this manually:
# python calc_inflowoutflow_stats.py -sf ../../02_precursor_shell/neutral_8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m/post_processing/abl_statistics129921.nc -ts 25000 -te 26800 -if setup_turbine_neutral.startAt25000.i

#rm -rf post_processing 
srun -n $cores --cpu_bind=cores $amrbin $input 2>&1

echo "Ending AMR-wind job at: " $(date)
