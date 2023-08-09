##################################
##################################
# AMR-Wind Turotial
##################################
##################################


##################################
##################################
# Compilation
##################################
##################################
git clone --recursive https://github.com/exawind/amr-wind.git
cd amr-wind/
mkdir build
cd build
 
cp ../../build.sh .
bash build.sh

# Now amr-wind has been compiled
# Look around the build directory and see the different files
# Let's run a test to make sure that compilation worked
# Source the environment first
source /nopt/nrel/ecom/exawind/exawind/scripts/amr-wind-gcc.sh
# Run all the unit tests
./amr_wind_unit_tests
# Run regression tests
ctest --output-on-failure -R act*

##################################
##################################
#  ABL
##################################
##################################
cd ABL
# Check the input file

# Request an allocation
salloc --time=60 --account=wks --ntasks=64

# Submit the simulation (should take 3-5 minutes to run)
bash runfile.sh

# On the browser go to:
# ed1.hpc.nrel.gov
# Open terminal on viz node and run:
source /nopt/nrel/ecom/exawind/exawind/scripts/amr-wind-gcc.sh
module load paraview/5.8.1-gui 
vglrun paraview


##################################
##################################
#  Wind turbine ALM
##################################
##################################
cd ALM

# Request an allocation
salloc --time=60 --account=wks --ntasks=64

# Submit the simulation (should take 10-15 minutes to run)
bash runfile.sh

# On the browser go to:
# ed1.hpc.nrel.gov
# Open terminal on viz node and run:
source /nopt/nrel/ecom/exawind/exawind/scripts/amr-wind-gcc.sh
module load paraview/5.8.1-gui
vglrun paraview


##################################
##################################
#  ABL with turbine
##################################
##################################
cd ABL_turbine

# Submit the simulation (should take 10-15 minutes to run)
bash runfile.sh

# On the browser go to:
# ed1.hpc.nrel.gov
# Open terminal on viz node and run:
source /nopt/nrel/ecom/exawind/exawind/scripts/amr-wind-gcc.sh
module load paraview/5.8.1-gui
vglrun paraview


