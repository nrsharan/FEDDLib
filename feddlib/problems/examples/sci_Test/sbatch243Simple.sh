#! /bin/bash -l

#SBATCH -N 4
#SBATCH --ntasks=216
#SBATCH -t 00:20:00
#SBATCH --output=216_Mono.out
#SBATCH --error=216_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_Test.exe --precfile=SIMPLE/parametersPrec_GDSW.xml
srun ./problems_sci_Test.exe --precfile=SIMPLE/parametersPrec_RGDSW.xml

