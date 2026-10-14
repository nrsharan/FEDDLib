# Artery drug transport (amlodipine study)

Structure-chemical interaction runs on five artery geometries: the artery is
pressurised, its smooth muscle activates, reorients and grows, and in the drug
runs from t = 860 s a drug enters through both walls — inner flags 5, 7, 8 and
outer flags 4, 6, 9 — and diffuses through the tissue. The control runs keep
the drug out (`Inflow Start Time` past the end of the run).

Build target `paper_amlodipine`. It needs Interface2 (the AceGen elements);
without it `main` returns `EXIT_FAILURE` after the MPI session, like every
other AceGen example here. The element is
`SMC_Interface_Active_Growth_Reorientation`; FEDDLib takes it with 34 history
values per Gauss point or with 35 (Interface2's `paper_amlodipine` branch, whose
dilution term keeps the previous step's det F as the 35th, starting at 1). The
runs take hours to days on a node and are therefore built but not registered as
ctests.

Every input file is a command-line option (`--simulationsParameters`,
`--materialParameters`, `--solverParameters`,
`--preconditionerParametersStructure`, `--preconditionerParametersChemistry`),
so all runs share the one executable.

## Inputs

- `base/`: the inputs the published results were computed with, one
  simulation and one material file per artery (dan, kim_guzman, narula,
  phinikaridou, plasschaert). The build directory gets dan's as
  `simulationParameters.xml` and `materialParameters_dan.xml`.
- `solverParameters.xml`, `preconditionerParameters_*.xml`: shared by all runs.
  The Newton iteration's linear solves use Eisenstat-Walker forcing terms
  (`Forcing Term Method` Type 2, tolerance between 1e-6 and 1e-3; the paper:
  Constant, 1e-5), which takes 15-17 % less time with results within about 1e-5
  (dan from 540 and 615 s, Fritz, 2026-10-12); the Newton tolerance is unchanged.
  Runs submitted before 2026-10-12 16:30 used Constant.
- `cases/`: one folder per run (90 in all) with its simulation and material
  file and `case.env` (`JOB_NAME`, `MATERIAL`, `MESH`, `NODES`, and `HOLD` for
  runs on hold). `cases/index.tsv` lists every run with its mesh, nodes, final
  time and what it changes. The folders are written by `make_cases.py` from
  `base/`; change the script, not the folders.
- The meshes are in `meshes/SPP2311` of the repository: the paper's
  `Artery_*_SCI.mesh` per artery (narula and phini: the `..._medium_SCI` ones)
  and the remeshed `narula_*`/`phini_*` of the mesh convergence study.

`make_cases.py` changes the paper's inputs in these ways for every run:

- **Adaptive time stepping.** The load ramp (0-1 s) keeps its fixed load steps.
  Every later time segment starts with the paper's `dt`, which is also its
  `Minimum dt`, and may grow to half the segment's length (`Maximum dt`); it
  doubles after 7 converged time steps (`Converged Time Steps Before Increase`;
  FEDDLib's default is 5).
  A time step in which an element cannot compute its state (its local Newton
  iteration fails) is repeated with a smaller time step. Only in the growth
  segment (540-840 s) is such a step accepted at the minimum time step, provided
  the global Newton iteration converges; anywhere else a failure at the minimum
  stops the run.
  The paper's runs (fixed time steps, an older FEDDLib) never stopped: every
  time step was accepted after at most 20 Newton iterations, converged or not,
  with no divergence test. A segment with `Accept Unconverged Newton At Minimum
  dt` set to true does that as a last resort: a time step that still fails at the
  minimum (also with element failures accepted) is repeated once more at that
  size without the divergence test and accepted after at most `MaxNonLinIts`
  iterations, unless the residual is not a finite number, ends above the residual
  the attempt started from, or the solver breaks down; the run then stops. (Without
  the divergence test a Newton iteration can blow up and then stop on a tiny
  update: the narula r0.1_L1 restart accepted ||F|| = 1.3e17 at 599.2 s with
  e03f42a and hung on the next step.) Each such step is reported in the `.out` file (`[adaptive] ... accepted
  the time step ... although the Newton iteration did not converge (||F|| = ...)`)
  and counted at the end. `make_cases.py` does not set it; it is switched on in a
  restart's `simulationParameters.xml` (growth segment) when a run cannot pass a
  time otherwise.
- **Output by time.** The solution and the postprocessing fields are written at
  the ends of the phases (1, 20, 220, 240, 540, 840, 860, 1000, 1200, 1500 s, and
  1600, 1700 s with unloading) and every 50 s (`Export Times`, `Export Interval`
  in `Exporter`). A time step ends exactly at the phase ends, and at the 50 s marks
  if it is not longer than 50 s; a longer one is written at its end.
- **Checkpoints** at 220, 540, 860 and 1200 s and every 25 s of growth (565 ... 840 s),
  in `checkpoints/`, so a run can be repeated from shortly before a time of interest.
- **The full stress tensor** (`Sxx` ... `Szz`) and the fibre stresses (`SFiber1`,
  `SFiber2`) among the postprocessing fields.
- **Preconditioner reuse.** The preconditioner is built in the first two Newton
  iterations of each time step and reused in the later ones (`Rebuild
  Preconditioner every Newton Iteration` false, `newtonLimit` 2 in `Parameter`;
  the paper: rebuilt in every iteration). With Eisenstat-Walker on, that takes
  10-14 % less time with results within about 1e-5 (dan from 540 and 615 s,
  faculty, 2026-10-12). It needs FEDDLib with the Newton iterations counted per
  time step (`NonLinearSolver`, 2026-10-12): an older one never rebuilds the
  preconditioner after the first time steps. Runs submitted before 2026-10-12
  21:00 rebuild it in every iteration.
- **phini**: `Alpha2` as in dan, kim and plasschaert (it was 10.0 in every region),
  and load steps of 0.01 s (with 0.02 s its first Newton iteration diverges).
  narula keeps the paper's `Alpha2` of 500 in its calcifications (19, 20), the
  lipid value, where the other arteries use 151.73775: as a volumetric penalty
  the higher value changes little.

The studies:

| Folder | Runs | Changes |
|---|---|---|
| `kappa_variation` | dan with three `Kappa` sets (15, 30, 45), control and drug | `Kappa` of media / degenerated media; 30 is the baseline |
| `arteries` | kim, narula, phini, plasschaert, control and drug | the artery |
| `pressure_drop_variation` | five arteries, drug, in `10_pressure_drop`, `20_pressure_drop`, `30_pressure_drop` | `Pressure Reduction Amount mmHg` 10, 20, 30 over 1000-1200 s |
| `diffusion_variation` | dan, drug | `D0` of adventitia / degenerated media (e = baseline) |
| `reaction_variation` | dan, drug | `M` of media / degenerated media (negative, as in the paper's files) |
| `degen_kappa_variation` | dan, control and drug | `Kappa` of the degenerated media at 33, 10, 0 % of the media's |
| `dose_variation` | dan, drug | `Inflow Concentration` 0.06-1.5 µM (paper: 2) |
| `drug_response_variation` | dan, drug | `C50` or `P` halved or doubled |
| `mesh_convergence` | narula and phini remeshed, control and drug | the mesh: the paper's geometry (`original`) and its cap shoulders smoothed with a 0.1 mm fillet (`r0.1`), each refined at the shoulders to 0.12, 0.06, 0.03 mm (L0-L2) |
| `axial_stretch` | dan, control and drug | `Axial Stretch` 0.05, ramped with the pressure |
| `basal_tone_tuning` | dan, control to 540 s | `Eta` × 0.75, 0.5, 0.35, 0.25, 0.15, 0.1 (all regions), to find the basal-tone levels |

The narula and phini runs on the paper's meshes are on hold (`HOLD` in
`case.env`) until `mesh_convergence` shows whether those meshes are fine
enough; `submit_cases.sh` skips them unless `INCLUDE_HELD=1` is set. The
larger meshes of `mesh_convergence` run on 2 or 4 nodes (`NODES` in
`case.env`).

The controls of the baseline arteries (`kappa_variation/artery_dan_30` and the
four in `arteries`) continue past 1500 s: the pressure goes to 0 over
1500-1600 s and is held to 1700 s.

Two parameters of `main.cpp` are new and default to the paper's behaviour:
`Inflow Concentration` (2.0) and `Axial Stretch` (0.0; the end faces and rings
get the axial displacement `Axial Stretch` × z, ramped up like the pressure).

## Output of a run

- `d_s`, `c` and `PostprocessingData` (`.h5` with `.xmf` for ParaView):
  displacement, concentration and the postprocessing fields at the output
  times.
- `time.txt` (the time of every accepted step), `newtonIterations.txt`,
  `linearIterations.txt`.
- `checkpoints/`: `Solution*.h5` (all checkpoint times) and `History<t>.h5`
  per checkpoint time.
- `ErrorLog.txt`: the messages of elements that could not compute their state,
  with element, Gauss point, time and time step. The element number is the
  0-based tetrahedron of the mesh file. Since Interface2 16ad68c the elements
  write to hidden per-rank files `.ErrorLog.raw.<pid>.txt` that are folded into
  `ErrorLog.txt` after every attempted time step and at the end of the run;
  FEDDLib does the fold when Interface2 has `newErrorLogStep()` (configure
  check `FEDD_HAVE_ACEGEN_ERRORLOG_STEP`). Most messages come from attempts that
  were repeated with a smaller time step; an attempt counts only if its `Time`
  and `DeltaT` match a step in `time.txt`.
- `<job>-<id>.out` / `.err`: the solver log, with a line `[adaptive] ...` for
  every repeated or accepted-with-failures time step.

## Running on Elysium

```
./submit_cases.sh [--dry-run] <build>/.../problems_paper_amlodipine.exe <run root> <case or study> ...
```

Each case runs in `<run root>/<case>`, which gets the case's input files, the
shared solver and preconditioner files, the mesh and `job.sh`; a case whose
directory already has a `job.sh` is skipped. `--dry-run` only prepares the
directories. `job.sh` asks for `NODES` exclusive nodes with 48 ranks each for
up to 7 days, loads `intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl`
(the executable must be built with the same modules), keeps every element
message in `ErrorLog.txt` (`INTERFACE2_ERRORLOG_EXAMPLES`; Interface2 otherwise
keeps 3 of each kind) and ends the job when one rank aborts
(`srun --kill-on-bad-exit=1`). `MAIL_USER=<address>` adds Slurm mails at the
end of each job; `INCLUDE_HELD=1` submits the runs on hold too.

## Continuing or repeating a run from a checkpoint

A restart writes its output files anew, so it runs in a new directory and
reads the checkpoints of the old one in place: copy the old run's input files,
mesh and `job.sh` to the new directory, and in its `simulationParameters.xml`
(`Timestepping Parameter`) set `Restart` to true, `Restart directory` to the
old run's `checkpoints` (e.g. `../artery_dan_30/checkpoints`), `Time step` to
the checkpoint time and `Import history` to true. Then submit `job.sh` there.
In the old run's directory the restart would overwrite `d_s`, `c`,
`PostprocessingData` and the text files from the start.
