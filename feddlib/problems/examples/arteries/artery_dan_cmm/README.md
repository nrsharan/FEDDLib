# artery_dan_cmm

The artery `meshes/SPP2311/Artery_dan_SCI.mesh` with the structure-chemistry
interaction (SCI) problem and the constrained-mixture (CMM) element, set up
as the FEDDLib side of svMultiPhysics's `tests/cases/def_diffu/artery_dan`:
both solve the same problem, so that the two codes and their FROSch
preconditioners can be compared one to one.

- Mesh: P1 mesh with seven tissue regions (volume flags 15-21); FEDDLib builds
  the P2 mesh (44300 nodes). `Artery_dan_SCI_svmp_flags.txt`, read through
  "Geometry Override" in the "Mesh Partitioner" list of
  `simulationParameters.xml`, gives the Dirichlet flags
  of exactly the nodes svMultiPhysics constrains (the flags are listed at the
  top of `main.cpp`).
- Material: every region uses svMultiPhysics's hollow_cylinder CMM parameter
  set (`materialParameters.xml`; placeholder parameters).
- Load: pressure on the inner wall ramped linearly to 85 mmHg at t = 1;
  zero concentration on the walls ("Inflow Start Time" 1e7).
- Time stepping: 50 steps of 0.02; Newmark for the displacement, backward
  Euler for the concentration.
- Linear solver: monolithic Belos Block GMRES (tolerance 1e-5, right
  preconditioning) with FROSch's `TwoLevelBlockPreconditioner`
  (`preconditionerParameters_Structure.xml`): RGDSW coarse space with the
  rotations for the displacement block, coarse basis entries below 1e-5
  dropped, the coarse basis of the first matrix kept ("Recycling",
  "Reuse: Coarse Basis"), no reuse of the symbolic factorization of the
  subdomain matrices, KLU2 everywhere. These are the fastest settings measured
  (below).

## Build and run on Elysium (HPC@RUB)

```
sbatch sampleConfigureScripts/elysium-build-job.sh problems_artery_dan_cmm
sbatch --ntasks=16 sampleConfigureScripts/elysium-run-example-job.sh \
    <build dir>/feddlib/problems/examples/arteries/artery_dan_cmm problems_artery_dan_cmm.exe
```

## Results (Elysium, cpu nodes, 50 steps)

All runs converge in every step. GMRES iterations are the mean per linear
solve; wall times of the same run vary by up to about 10% between nodes.

| Run | 16 ranks: wall, GMRES its, Newton its/step | 32 ranks |
|---|---|---|
| FEDDLib, these settings | 3233 s, 203.7, 3.8 | 1735 s, 266.0, 3.8 |
| FEDDLib, coarse basis recomputed for every matrix | 3464 s, 156.5, 3.8 | 1545 s, 192.2, 3.8 |
| svMultiPhysics `trilinos-frosch-block`, the same FROSch settings | 1608 s, 26.8, 3 | 702 s, 34.7, 3 |

- Recomputing the coarse basis lowers the iterations but takes about as long
  overall (7% slower on 16 ranks, 11% faster on 32).
- Without the rotations, FEDDLib needs slightly more iterations (first time
  step, 16 ranks: 199, 164, 181, 163 against 171, 146, 160, 142).
- Where the time goes (32 ranks, coarse basis recomputed): FROSch `compute`
  724 s for 191 matrices, GMRES 577 s. FROSch `initialize` (once) takes 0.4 s.

## Comparison with svMultiPhysics

The first Newton system of both codes was written (`FEDD_WRITE_SYSTEM`, below;
`SVMP_TRILINOS_WRITE_SYSTEM` in svMultiPhysics) and compared with
svMultiPhysics's `tests/cases/def_diffu/tools/compare_systems.py`:

- The displacement blocks agree up to a constant factor (0.1%); the
  concentration blocks differ by the time integration (30%).
- The concentration is zero, so the displacement-concentration blocks are
  zero: FEDDLib stores none, svMultiPhysics stored them as zeros.
- FEDDLib replaces the Dirichlet rows by rows of the identity and keeps the
  Dirichlet columns; svMultiPhysics removes both and scales the system.

svMultiPhysics needs about 6 times fewer GMRES iterations on the same
problem, also with FEDDLib's FROSch settings. The difference is already in
FROSch's first level: with `OneLevelPreconditioner` ("FROSch Preconditioner
Type"), FEDDLib needs 557, 477, 530, 470 iterations in the first time step on
16 ranks and svMultiPhysics 55, 52, 52. Not the cause: the systems, the
Dirichlet treatment and scaling (emulated in svMultiPhysics), the FROSch
settings, the rotations, the Jacobian given to GMRES (it equals the system
FROSch is built from in every Newton iteration) and the partitions (connected,
of about the same size). Left (svMultiPhysics) against right (FEDDLib)
preconditioning and how the two codes hand the matrix and maps to FROSch
remain untested.

## Debugging output (FEDD_WRITE_SYSTEM)

With the environment variable `FEDD_WRITE_SYSTEM=<prefix>`:

- the monolithic preconditioner writes the blocks of the first system it is
  built for (`<prefix>ij.mm`, MatrixMarket; `Preconditioner_def.hpp`);
- this example writes the P2 nodes of each process
  (`<prefix>_coords_<rank>.txt`: node GID and coordinates);
- `TimeProblem` compares NOX's Jacobian with the merged system it is copied
  from in every Newton iteration and prints the result ("W check").
