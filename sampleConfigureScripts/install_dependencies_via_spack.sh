#!/bin/bash

# Install Dependencies via spack

spack install hdf5%oneapi build_type=RelWithDebInfo ^intel-oneapi-mpi

spack install boost %oneapi +mpi cxxstd=17 cppflags='-lpthread' cxxflags='-lpthread' cflags='-lpthread'  ^intel-oneapi-mpi

spack install metis %oneapi build_type=RelWithDebInfo +int64 +real64

# Note that the last part is the hash of the installed metis library. The hash can be obtained using `spack find -l`
spack install parmetis %oneapi +int64 build_type=RelWithDebInfo ^intel-oneapi-mpi ^/vyzvexk
