#!/bin/bash -l
#SBATCH --job-name=feddlib-build
#SBATCH --account=balzadlb_0000
#SBATCH --partition=cpu_filler
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/lustre/nurans63/feddlib-logs/%x-%j.out
#
# Configures (if needed) and builds FEDDLib on Elysium (HPC@RUB) on a compute
# node (computations, compiling included, are not allowed on the login
# nodes), with the Intel oneAPI compilers, the Trilinos in
# ~/opt/intel-llvm/trilinos (global ordinal long long, with FROSch and
# PARDISO_MKL) and an Interface2 install with the constrained-mixture (CMM)
# element. The flags are those of the earlier Intel build
# (/lustre/nurans63/feddlib_comparison/do-config-feddlib_intel.sh).
#
#   mkdir -p /lustre/nurans63/feddlib-logs
#   sbatch sampleConfigureScripts/elysium-build-job.sh [make targets, e.g. problems_artery_dan_cmm]
# (TriBITS names an example's target <package>_<example>; its executable is
# <package>_<example>.exe in the example's build directory.)
#
# The paths are set below; pass the make targets as arguments rather than
# variables with sbatch --export=ALL,... (that copies the submitting shell's
# environment, whose module setup may be incomplete, into the job).
set -eo pipefail
unset SLURM_EXPORT_ENV

SOURCE_DIR=$HOME/dev/FEDDLib/FEDDLib-cmm
BUILD_DIR=/lustre/nurans63/feddlib_artery_cmm
TRILINOS_DIR=$HOME/opt/intel-llvm/trilinos
INTERFACE_DIR=$HOME/opt/intel-llvm/interface2-latest
BUILD_TYPE=RelWithDebInfo

ml load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl cmake

FLAGS="-D HAVE_EXPLICIT_INSTANTIATION -D FROSCH_Epetra64 -D WeUseTpetra -Wno-deprecated -Wno-sign-compare -Wno-unused-variable  -D NEW_PARMETIS -D NEW_METIS -fpermissive -D FEDD_HAVE_ACEGENINTERFACE"

echo "== $(date)  $(hostname)  source $SOURCE_DIR  build $BUILD_DIR"
mkdir -p $BUILD_DIR
cd $BUILD_DIR
if [ ! -f CMakeCache.txt ]; then
  cmake \
  -G "Unix Makefiles" \
  -D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
  -D CMAKE_C_COMPILER=mpiicx \
  -D CMAKE_CXX_COMPILER=mpiicpx \
  -D CMAKE_Fortran_COMPILER=mpiifx \
  -D CMAKE_C_FLAGS:STRING="-fiopenmp $FLAGS" \
  -D CMAKE_CXX_FLAGS:STRING="-fp-model=precise $FLAGS -Wno-write-strings -Wno-cpp -w" \
  -D CMAKE_Fortran_FLAGS:STRING="-fiopenmp" \
  -D CMAKE_CXX_STANDARD:STRING=14 \
  -D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
  -D FEDDlib_ENABLE_ALL_PACKAGES:BOOL=ON \
  -D FEDDlib_ENABLE_TESTS:BOOL=ON \
  -D TPL_FIND_SHARED_LIBS:BOOL=OFF \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D TPL_ENABLE_Trilinos:BOOL=ON \
  -D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
  -D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib64 \
  -D TPL_ENABLE_AceGENInterface:BOOL=ON \
  -D TPL_AceGENInterface_LIBRARIES:STRING="$INTERFACE_DIR/lib64/libinterface2.a;" \
  -D TPL_AceGENInterface_INCLUDE_DIRS:STRING="$INTERFACE_DIR/include" \
  $SOURCE_DIR
fi
make -j ${SLURM_CPUS_PER_TASK:-8} "$@"
echo "== done $(date)"
