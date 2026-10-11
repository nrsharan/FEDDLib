# Artery drug transport (amlodipine study)

Structure-chemical interaction runs on five artery geometries: the artery is
pressurised, its smooth muscle activates, reorients and grows, and from
t = 860 s a drug enters through both walls — inner flags 5, 7, 8 and outer
flags 4, 6, 9 — and diffuses through the tissue.

Build target `paper_amlodipine`. It needs Interface2 (the AceGen elements);
without it `main` returns `EXIT_FAILURE` after the MPI session, like every
other AceGen example here. The runs take days on a node and are therefore
built but not registered as ctests.

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
- `cases/`: one folder per run with its simulation and material file and
  `case.env` (job name, material file, mesh). `cases/index.tsv` lists every run
  and what it changes. The folders are written by `make_cases.py` from `base/`;
  change the script, not the folders.

`make_cases.py` changes the paper's inputs in these ways for every run:

- **Adaptive time stepping.** The load ramp (0-1 s) keeps its fixed load steps.
  Every later time segment starts with the paper's `dt`, which is also its
  `Minimum dt`, and may grow to half the segment's length (`Maximum dt`).
  Elements whose local Newton iteration fails are accepted at the minimum time
  step in the growth segment (540-840 s) only; anywhere else such a step is
  repeated with a smaller time step, and a failure at the minimum stops the run.
- **Output by time.** The solution and the postprocessing fields are written at
  the ends of the phases (1, 20, 220, 240, 540, 840, 860, 1000, 1200, 1500 s, and
  1600, 1700 s with unloading) and every 50 s (`Export Times`, `Export Interval`
  in `Exporter`). A time step ends exactly at the phase ends, and at the 50 s marks
  if it is not longer than 50 s; a longer one is written at its end.
- **Checkpoints** at 220, 540, 860 and 1200 s and every 25 s of growth (565 ... 840 s),
  in `checkpoints/`, so a run can be repeated from shortly before a time of interest.
- **The full stress tensor** (`Sxx` ... `Szz`) and the fibre stresses (`SFiber1`,
  `SFiber2`) among the postprocessing fields.
- **phini**: `Alpha2` as in dan, kim and plasschaert (it was 10.0 in every region),
  and load steps of 0.01 s (with 0.02 s its first Newton iteration diverges).

The studies:

| Folder | Runs | Changes |
|---|---|---|
| `kappa_variation` | dan with three `Kappa` sets (15, 30, 45), control and drug | `Kappa` of media / degenerated media; 30 is the baseline |
| `arteries` | kim, narula, phini, plasschaert, control and drug | the artery |
| `pressure_drop_variation` | five arteries, drug | `Pressure Reduction Amount mmHg` 10, 20, 30 over 1000-1200 s |
| `diffusion_variation` | dan, drug | `D0` of adventitia / degenerated media (e = baseline) |
| `reaction_variation` | dan, drug | `M` of media / degenerated media |
| `degen_kappa_variation` | dan, control and drug | `Kappa` of the degenerated media at 33, 10, 0 % of the media's |
| `dose_variation` | dan, drug | `Inflow Concentration` 0.06-1.5 µM (paper: 2) |
| `drug_response_variation` | dan, drug | `C50` or `P` halved or doubled |
| `mesh_convergence` | narula and phini remeshed, control and drug | the mesh: the paper's geometry (`original`) and its cap shoulders smoothed with a 0.1 mm fillet (`r0.1`), each refined at the shoulders to 0.12, 0.06, 0.03 mm (L0-L2) |
| `axial_stretch` | dan, control and drug | `Axial Stretch` 0.05, ramped with the pressure |
| `basal_tone_tuning` | dan, control to 540 s | `Eta` scaled, to find the basal-tone levels |

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

## Running on Elysium

```
./submit_cases.sh <build>/.../problems_paper_amlodipine.exe <run root> <case or study> ...
```

Each case runs in `<run root>/<case>` with 48 ranks per node. `--dry-run`
only prepares the directories.

A run stopped by the time limit continues from its last checkpoint: in its
`simulationParameters.xml` set `Restart` to true, `Restart directory` to
`checkpoints`, `Time step` to the checkpoint time and `Import history` to true,
and submit `job.sh` again.
