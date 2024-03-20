#module purge
#module load comp-intel mkl cmake
#export CC=icc
#export FC=ifort
#export CXX=icpc
#CC=icc
#FC=ifort
#CXX=icpc

module purge
module load craype-x86-spr
#module load intel-oneapi-mkl/2023.2.0-intel
#module load intel-oneapi-mpi/2021.10.0-intel
#module load intel-oneapi-compilers/2023.2.0
module load intel-oneapi-mkl
module load intel-oneapi-compilers
module load cmake git



