#!/bin/bash -l

module purge
module load gcc
module load mpt
module load cmake
#module load netcdf-c/4.7.3

export EXAWIND_DIR=/nopt/nrel/ecom/exawind/exawind-2020-09-21/install/gcc

cmake -DAMR_WIND_ENABLE_CUDA:BOOL=OFF \
      -DAMR_WIND_ENABLE_MPI:BOOL=ON \
      -DAMR_WIND_ENABLE_OPENMP:BOOL=OFF \
      -DAMR_WIND_TEST_WITH_FCOMPARE:BOOL=OFF \
      -DCMAKE_BUILD_TYPE=Release \
      -DAMR_WIND_ENABLE_NETCDF:BOOL=ON \
      -DNETCDF_DIR:PATH=/nopt/nrel/ecom/hpacf/software/2020-07/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.4.0/netcdf-c-4.7.3-533s5vfhvbbvpgxambbzk66vtlcce2u6  \
      -DAMR_WIND_ENABLE_OPENFAST:BOOL=ON \
      -DOpenFAST_ROOT:PATH=${EXAWIND_DIR}/openfast \
      -DAMR_WIND_ENABLE_HYPRE:BOOL=OFF \
      -DHYPRE_ROOT:PATH=${EXAWIND_DIR}/hypre \
      -DAMR_WIND_ENABLE_MASA:BOOL=OFF \
      -DAMR_WIND_ENABLE_TESTS:BOOL=ON \
      -DAMR_WIND_ENABLE_FORTRAN:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -DAMR_WIND_ENABLE_ALL_WARNINGS:BOOL=ON \
      -DBUILD_SHARED_LIBS:BOOL=ON \
      ..

nice make -j16
