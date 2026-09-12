#!/bin/bash
# FEDDLib (this fork's CMake: project FEDDlib) on macOS / Apple silicon (M4)
# against the Trilinos 16.1.0 from configure-Trilinos-on-macOS-AppleSilicon.sh
# (Homebrew Open MPI, Apple Accelerate via veclibfort, Open MPI-built spack
# HDF5/ParMETIS), with the Interface2 AceGen interface, as used for the local
# builds of 2026-09. The FEDDLib develop line uses other CMake names (project
# FEDDLib, HDF5 and zlib as TPLs); see its doc/sampleConfigureScripts.
#
# Build order: install-dependencies-on-macOS-AppleSilicon-via-spack.sh,
# configure-Trilinos-on-macOS-AppleSilicon.sh, Interface2 (its
# config_scripts/do-config-interface2-macos-accelerate.sh, installed to
# ~/opt/interface2), then this script.
# Usage: run from an empty build directory, then `ninja`. Trilinos installs
# its libraries in lib/ here (lib64/ on the cluster). Without MKL, replace
# "pardisomkl" by "klu2" in the parameter files.

BUILD_TYPE=${BUILD_TYPE:-RelWithDebInfo}
#BUILD_TYPE=Debug

SOURCE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TRILINOS_DIR=${TRILINOS_DIR:-$HOME/opt/trilinos-16.1.0}
INSTALL_DIR=${INSTALL_DIR:-$HOME/opt/FEDDLib}
INTERFACE_DIR=${INTERFACE_DIR:-$HOME/opt/interface2}

rm -rf CMake*
cmake \
-G Ninja \
-D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
-D CMAKE_C_COMPILER=mpicc \
-D CMAKE_CXX_COMPILER=mpicxx \
-D CMAKE_Fortran_COMPILER=mpifort \
-D CMAKE_INSTALL_PREFIX:STRING=${INSTALL_DIR} \
-D CMAKE_C_FLAGS:STRING="-D HAVE_EXPLICIT_INSTANTIATION -D FROSCH_Epetra64 -D WeUseTpetra -Wno-deprecated -Wno-sign-compare -Wno-unused-variable  -D NEW_PARMETIS -D NEW_METIS -fpermissive -D FEDD_HAVE_ACEGENINTERFACE" \
-D CMAKE_CXX_FLAGS:STRING="-D HAVE_EXPLICIT_INSTANTIATION -D FROSCH_Epetra64 -D WeUseTpetra -Wno-deprecated -Wno-sign-compare -Wno-unused-variable  -D NEW_PARMETIS -D NEW_METIS -fpermissive -D FEDD_HAVE_ACEGENINTERFACE -Wno-write-strings -Wno-cpp -w" \
-D CMAKE_Fortran_FLAGS:STRING="" \
-D CMAKE_CXX_STANDARD:STRING=14 \
-D MPI_EXEC_MAX_NUMPROCS:STRING=6 \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D FEDDlib_ENABLE_ALL_PACKAGES:BOOL=ON \
-D FEDDlib_ENABLE_TESTS:BOOL=ON \
-D TPL_FIND_SHARED_LIBS:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D TPL_ENABLE_Trilinos:BOOL=ON \
-D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
-D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib \
-D TPL_ENABLE_AceGENInterface:BOOL=ON \
-D TPL_AceGENInterface_LIBRARIES:STRING="$INTERFACE_DIR/lib/libinterface2.a;" \
-D TPL_AceGENInterface_INCLUDE_DIRS:STRING="$INTERFACE_DIR/include" \
${SOURCE_DIR}
