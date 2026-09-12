#!/bin/bash

BUILD_TYPE=DEBUG

TRILINOS_DIR=`spack location -i trilinos%gcc`
SOURCE_DIR=$HOME/dev/FEDDLib/FEDDLib-Sharan/FEDDLib
INSTALL_DIR=$FASTTMP/opt/intel-llvm/FEDDLib
INTERFACE_DIR=$HOME/opt/gcc/Interface2-develop-debug
# -D CMAKE_CXX_STANDARD_LIBRARIES:STRING="/opt/local/lib/gcc9/libgfortran.dylib /opt/local/lib/mpich-gcc9/libmpifort.dylib" \

rm -rf CMake*
cmake \
-G Ninja \
-D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
-D CMAKE_C_COMPILER=/home/sramesh/opt/spack/opt/spack/linux-manjaro24-skylake/gcc-14.1.1/openmpi-5.0.3-ei3fomn2pktkdjvzgbgx5wd4ckswmngt/bin/mpicc \
-D CMAKE_CXX_COMPILER=/home/sramesh/opt/spack/opt/spack/linux-manjaro24-skylake/gcc-14.1.1/openmpi-5.0.3-ei3fomn2pktkdjvzgbgx5wd4ckswmngt/bin/mpic++ \
-D CMAKE_Fortran_COMPILER=/home/sramesh/opt/spack/opt/spack/linux-manjaro24-skylake/gcc-14.1.1/openmpi-5.0.3-ei3fomn2pktkdjvzgbgx5wd4ckswmngt/bin/mpifort \
-D CMAKE_INSTALL_PREFIX:STRING=${INSTALL_DIR} \
-D CMAKE_C_FLAGS:STRING="-fopenmp -D HAVE_EXPLICIT_INSTANTIATION -D FROSCH_Epetra64 -D WeUseTpetra -Wno-deprecated -Wno-sign-compare -Wno-unused-variable  -D NEW_PARMETIS -D NEW_METIS -fpermissive -D FEDD_HAVE_ACEGENINTERFACE" \
-D CMAKE_CXX_FLAGS:STRING="-D HAVE_EXPLICIT_INSTANTIATION -D FROSCH_Epetra64 -D WeUseTpetra -Wno-deprecated -Wno-sign-compare -Wno-unused-variable  -D NEW_PARMETIS -D NEW_METIS -fpermissive -D FEDD_HAVE_ACEGENINTERFACE" \
-D CMAKE_Fortran_FLAGS:STRING="-fopenmp" \
-D CMAKE_CXX_STANDARD:STRING=14 \
-D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D FEDDlib_ENABLE_ALL_PACKAGES:BOOL=ON \
-D FEDDlib_ENABLE_TESTS:BOOL=ON \
-D TPL_FIND_SHARED_LIBS:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D TPL_ENABLE_Trilinos:BOOL=ON \
-D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
-D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib \
-D TPL_ENABLE_AceGENInterface:BOOL=ON \
-D TPL_AceGENInterface_LIBRARIES:STRING="$INTERFACE_DIR/lib/libinterface2.a;$INTERFACE_DIR/lib/libaceutility.a;" \
-D TPL_AceGENInterface_INCLUDE_DIRS:STRING="$INTERFACE_DIR/include" \
${SOURCE_DIR}
