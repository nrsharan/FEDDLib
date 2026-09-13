# ccb_hollow_cylinder

A hollow cylinder (Ri=1.0, Ro=1.3, L=2.0) with the structure-chemistry
interaction (SCI) problem and the non-CMM SMC element
(`SCI_SMC_Active_Growth_Reorientation`), solved with FEDDLib's
NOX/Belos/FROSch solver stack. `../ccb_hollow_cylinder_cmm` is the same case
with the constrained-mixture (CMM) element.

- Boundary conditions: bottom and top faces fixed axially, three
  circumferential pins for static determinacy, pressure on the inner wall
  ramped 0 -> 20 (plain pressure units) over t = [0, 10].
- All `*Bool` material flags are 0 (purely passive, load-driven response).
- `solverParameters.xml` prints NOX's outer iterations (`Outer Iteration`,
  `Outer Iteration StatusTest`) and Belos's details of every linear solve
  (`Linear Solver Details`, `Inner Iteration`).

## Build

From a configured Trilinos+FEDDLib(+Interface2/AceGenInterface) build
directory:

```
make ccb_hollow_cylinder
# or: make -j <N>   (builds everything, this target included)
```

## Run

```
cd <build_dir>/feddlib/problems/examples/arteries/ccb_hollow_cylinder
mpirun -np 4 ./ccb_hollow_cylinder.exe 2>&1 | tee run.log
```

(The binary may be called `problems_ccb_hollow_cylinder.exe` depending on
TriBITS target aliasing in this build -- check
`ctest -N | grep ccb_hollow_cylinder` if `ccb_hollow_cylinder.exe` isn't
found next to the copied XML/mesh files.)

## Mesh / BC flags

See `meshes/ccb_hollow_cylinder/hollow_cylinder_p1.mesh` (P1; FEDDLib builds
the P2 mesh) and the comment block at the top of `main.cpp` for the
face/vertex flags (bottom/top axial fix, 3 circumferential pins, inner-wall
pressure).
