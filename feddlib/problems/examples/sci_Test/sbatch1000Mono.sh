#! /bin/bash -l

#SBATCH -N 14
#SBATCH --ntasks=1000
#SBATCH -t 00:20:00
#SBATCH --output=1000_Mono.out
#SBATCH --error=1000_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_Test.exe --precfile=parametersPrec_GDSW.xml	--problemfile=parametersProblemSCI_GDSW.xml
srun ./problems_sci_Test.exe --precfile=parametersPrec_RGDSW.xml --problemfile=parametersProblemSCI_RGDSW.xml


