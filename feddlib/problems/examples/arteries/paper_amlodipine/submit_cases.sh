#!/bin/bash
# Prepares a run directory for each given case of cases/ and submits it to Slurm (Elysium).
#
#   submit_cases.sh [--dry-run] <executable> <run root> <case> [<case> ...]
#
# <executable>  the built problems_paper_amlodipine.exe
# <run root>    where the runs go; a case runs in <run root>/<case>
# <case>        a folder under cases/, e.g. kappa_variation/artery_dan_withDrug_30,
#               or a whole study, e.g. pressure_drop_variation
#
# The run directory gets the case's two input files, the shared solver and preconditioner
# files, the mesh, and job.sh. --dry-run prepares the directories without submitting.
# MAIL_USER=<address> in the environment adds Slurm mails at the end of each job.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
MESHES=$(cd "$HERE/../../../../../meshes/SPP2311" && pwd)

DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then DRY_RUN=1; shift; fi
if [ $# -lt 3 ]; then sed -n '2,15p' "$0"; exit 1; fi
EXE=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
RUN_ROOT=$2
shift 2
[ -x "$EXE" ] || { echo "Not an executable: $EXE"; exit 1; }

CASES=()
for arg in "$@"; do
    arg=${arg%/}
    [ -d "$HERE/cases/$arg" ] || { echo "No case or study cases/$arg"; exit 1; }
    while IFS= read -r env; do CASES+=("$(dirname "${env#$HERE/cases/}")"); done < <(find "$HERE/cases/$arg" -name case.env | sort)
done

for case in "${CASES[@]}"; do
    source "$HERE/cases/$case/case.env"
    run="$RUN_ROOT/$case"
    if [ -e "$run/job.sh" ]; then echo "skipped (exists): $run"; continue; fi
    mkdir -p "$run"
    cp "$HERE/cases/$case/simulationParameters.xml" "$HERE/cases/$case/$MATERIAL" "$run/"
    cp "$HERE/solverParameters.xml" "$HERE/preconditionerParameters_Structure.xml" "$HERE/preconditionerParameters_Chemistry.xml" "$run/"
    cp "$MESHES/$MESH" "$run/"
    mail=""
    [ -n "${MAIL_USER:-}" ] && mail="#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$MAIL_USER"
    cat > "$run/job.sh" <<EOF
#!/bin/bash -l
#SBATCH --job-name=$JOB_NAME
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-core=1
#SBATCH --exclusive
#SBATCH --account=balzadlb_0000
#SBATCH --partition=cpu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
$mail

unset SLURM_EXPORT_ENV
# the modules the executable was built with (build-job.sh)
ml load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl

srun --mpi=pmi2 $EXE --materialParameters=$MATERIAL
EOF
    if [ $DRY_RUN -eq 1 ]; then
        echo "prepared: $run"
    else
        (cd "$run" && echo "$case: $(sbatch job.sh)")
    fi
done
