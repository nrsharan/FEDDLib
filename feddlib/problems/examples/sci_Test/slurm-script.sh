#!/bin/bash -l
#
#SBATCH --job-name=partialArtery
#SBATCH --comment="Partial artery"
#SBATCH --time=12:00:00         ### time the job will appr. run
#SBATCH --nodes=3
#SBATCH --output=/lustre/k105be/k105be13/slurm/outputs/partialArtery.%j.out              ### output file for console output
#SBATCH --error=/lustre/k105be/k105be13/slurm/outputs/partialArtery.%j.err               ### output file for console error
#SBATCH --ntasks=144              ### Number of tasks per job    (or next line, should one of them)
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV

source ~/opt/spack/share/spack/setup-env.sh
#spack load intel-oneapi-compilers
spack load intel-oneapi-mkl
spack load hdf5
spack load metis
spack load parmetis
spack load boost
#spack load cmake
#export I_MPI_PMI_LIBRARY=/lib64/libpmi2.so
spack --version
#which mpiexec
srun ./problems_sci_Test.exe
#mpirun -np 216 ./problems_sci_Test.exe
