#!/bin/bash -l
#
#SBATCH --job-name=artery_kim_withDrug
#SBATCH --comment="Artery_kim_withDrug"
#SBATCH --time=24:00:00         ### time the job will appr. run
#SBATCH --nodes=3               ### Node count required for the job
#SBATCH --output=/home/hpc/k105be/k105be13/slurm/outputs/kim_artery_withDrug.%j.out              ### output file for console output
#SBATCH --error=/home/hpc/k105be/k105be13/slurm/outputs/kim_artery_withDrug.%j.err               ### output file for console error
#SBATCH --ntasks=108              ### Number of tasks per job    (or next line, should one of them)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharan.nuraniramesh@rub.de

unset SLURM_EXPORT_ENV

module load intel/2023.2.1 intelmpi/2021.10.0 mkl/2023.2.0
srun ./problems_artery_kim_withDrug.exe