#!/bin/bash -l
#SBATCH --job-name=feddlib-run
#SBATCH --account=balzadlb_0000
#SBATCH --partition=cpu_filler
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-core=1
#SBATCH --output=/lustre/nurans63/feddlib-logs/%x-%j.out
#
# Runs a built FEDDLib example on Elysium (HPC@RUB) in a scratch copy of its
# build directory under /lustre, then summarizes the Newton and GMRES
# iterations:
#   sbatch --ntasks=16 sampleConfigureScripts/elysium-run-example-job.sh <example build dir> <executable> [final time]
# e.g.
#   sbatch --ntasks=16 sampleConfigureScripts/elysium-run-example-job.sh \
#       /lustre/nurans63/feddlib_artery_cmm/feddlib/problems/examples/arteries/artery_dan_cmm problems_artery_dan_cmm.exe 0.2
# A final time replaces "Final time" in simulationParameters.xml. Pass options
# as arguments, not with sbatch --export=ALL,... (see elysium-build-job.sh).
set -o pipefail
unset SLURM_EXPORT_ENV

EXAMPLE_DIR=${1:?usage: $0 <example build dir> <executable> [final time]}
EXE=${2:?usage: $0 <example build dir> <executable> [final time]}
FINAL_TIME=${3:-}
RUN_ROOT=/lustre/nurans63/feddlib-runs
NP=${SLURM_NTASKS:-1}

ml load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl cmake

T=$RUN_ROOT/$(basename $EXAMPLE_DIR)_np${NP}_${SLURM_JOB_ID:-local}
mkdir -p $T
cd $EXAMPLE_DIR
cp *.xml *.txt *.mesh $T/ 2>/dev/null
ln -sf $EXAMPLE_DIR/$EXE $T/$EXE
cd $T
if [ -n "$FINAL_TIME" ]; then
  sed -i -E "s|(<Parameter name=\"Final time\" type=\"double\" value=\")[^\"]*\"|\1$FINAL_TIME\"|" simulationParameters.xml
fi
echo "=== $(basename $EXAMPLE_DIR)  np=$NP  final time=$(grep -o 'name="Final time" type="double" value="[^"]*"' simulationParameters.xml | grep -o '[0-9.e+-]*"$' | tr -d '"')  node=$(hostname)  dir=$T"

start=$(date +%s)
srun --mpi=pmi2 ./$EXE > run.log 2>&1
rc=$?
echo "EXIT $rc  wall $(( $(date +%s) - start )) s"
echo "errors/exceptions: $(grep -a -c -i 'error\|terminat\|exception' run.log)"
# FEDDLib reports, per time step, the number of Newton steps and the average
# number of GMRES iterations of their linear solves ("Number of Iterations" in
# the NOX status tests is the Newton count, not GMRES).
echo "time steps: $(grep -a -c 'Total nonlinear iterations' run.log)"
echo "Newton steps per time step (min/mean/max) and GMRES iterations per linear solve (min/mean/max of the time-step averages):"
grep -a "Total nonlinear iterations" run.log | sed -E 's/.*iterations : *([0-9]+) +with an average of +([0-9.]+) linear.*/\1 \2/' | \
  awk '{n++; s1+=$1; s2+=$2; if(a==""||$1<a)a=$1; if($1>b)b=$1; if(c==""||$2<c)c=$2; if($2>d)d=$2}
       END{if(n) printf "  %d / %.1f / %d    %.1f / %.1f / %.1f\n", a, s1/n, b, c, s2/n, d}'
grep -a "Total nonlinear iterations" run.log | sed 's/^/  /'
