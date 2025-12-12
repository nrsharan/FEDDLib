#!/bin/bash -l
#
#SBATCH --job-name=plass
#SBATCH --comment="Plass"
#SBATCH --time=7-00:00:00         ### time the job will appr. run
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --output=/home/nurans63/slurm/outputs/%x-%j.out              ### output file for console output
#SBATCH --error=/home/nurans63/slurm/outputs/%x-%j.err               ### output file for console error
#SBATCH --ntasks=48              ### Number of tasks per job    (or next line, should one of them)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharan.nuraniramesh@rub.de
#SBATCH --account=balzadlb_0000
#SBATCH --partition=cpu
#SBATCH --ntasks-per-core=1
#SBATCH --exclusive

unset SLURM_EXPORT_ENV

module load intel-oneapi-mkl
module load intel-oneapi-mpi/2021.12.1-qdmj2yh
srun --mpi=pmi2 ./problems_artery_plass_pd_0.exe