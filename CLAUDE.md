# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A CFD (Computational Fluid Dynamics) research codebase for an undergraduate thesis (TCC) on
wind loading simulation. It solves the incompressible Navier-Stokes equations with an explicit
Taylor-Galerkin finite-element scheme on 8-node hexahedral meshes, accelerated on GPU with
CUDA Fortran. Typical cases: flow over buildings/terrain, boundary-layer profiles, immersed
rigid bodies, wind-tunnel-style cavity flows (see `mshs/*.geo` and `artigos/*.pdf` for the
underlying literature — Braun, Madalozzo, Chorin pseudo-compressibility, Kim & Hangan, etc.).

There is no test suite — this is a numerical simulation / research pipeline, validated by
running cases and comparing output fields, not by unit tests.

## Pipeline / workflow

1. **Mesh**: build geometry + mesh in Gmsh (`mshs/*.geo` → `make mesh` → `.msh`, format msh2,
   version-2.2 ASCII — the Fortran/Python readers depend on this format).
2. **Preprocess**: turn the `.msh` + a `.phys` (boundary-condition spec) into the raw numeric
   files the solver reads. Two parallel implementations exist:
   - `src/TG2_preprocessor.py` — current, Python + `gmsh` API. Entry point:
     `python src/TG2_preprocessor.py <mesh.msh> <pname> [--show] [--dir DIR] [--config file.yaml]`.
     `pname` must be **exactly 6 characters, no spaces** — it's used as the filename prefix for
     every generated file. Without `--config` it interactively prompts (Portuguese) for boundary
     conditions/parameters/fluid properties per physical group and saves the answers to a
     timestamped `config_*.yaml` for reuse via `--config`.
   - `src/interface.f90` — legacy Fortran preprocessor, reads `<pname>.phys`, the `.msh`, and
     `INTERFACE.dat` (canned answers) instead of prompting. Build with `make interface`.
   - Both produce the same file set: `<pname>000.cor` (node coords), `<pname>.con`
     (hex connectivity), `<pname>CC.bv`/`CC.bp` (velocity/pressure BCs), `<pname>000.v`/`000.pr`
     (initial fields), `<pname>CS##.nnn`/`CS##.sup` (immersed-object normals/faces, one pair per
     body), `inicia.par` (run-control parameters), `<pname>.pro` (fluid properties + immersed
     object geometry), and `PARAMETER.dat` (a Fortran `PARAMETER(...)` block with mesh-size
     constants to paste into the solver).
3. **Solve**: `src/cuda_cfd.cuf`, program `CUDA_HEXA_IGOR`, built with `make cfd` (`nvfortran`,
   requires the NVIDIA HPC SDK and a CUDA GPU). Near the top of the program, a hardcoded
   `PARAMETER (NEMAX=..., NNOS=..., MNOBJ=..., NNCSX=..., NFCSX=..., NBU=..., NBV=..., NBW=...,
   NBP=..., NCMAX=..., NNMAX=...)` block **must be edited to match the current mesh/case**
   (copy the values from the `PARAMETER.dat` generated in the preprocessing step) before
   compiling — these are compile-time array dimensions, not runtime input. Old parameter sets
   for past cases are left commented out above the active one; follow that convention when
   swapping cases. `src/FLUexp2p.f90` is the original serial CPU reference implementation
   (by Alexandre Braun) that the CUDA port is adapted from — consult it to understand intended
   numerical behavior when the CUDA version is unclear, but it is not built by the Makefile.
4. **Run**: on the cluster via SLURM (`tupi.slurm`, partition `tupi`, 1 GPU) or
   `expr/*/submit.slurm` variants. The solver reads the preprocessed files from its working
   directory and writes timestepped `<pname>NNN.v` / `<pname>NNN.pr` (and `.cor`/`.con` once).
5. **Postprocess**: `src/make_vtk.f90` (`make vtk` → `./vtks <pname> <iter>`) converts one
   timestep's `.cor`/`.con`/`.v`/`.pr` into a legacy ASCII `.vtk` for ParaView. Batch-convert a
   range with `src/build_vtks.sh <pname> <first_iter> <num_iters>` (loops calling `./vtks`).
   `src/make_vtk_media.f90` is the equivalent for time-averaged ("media") fields.

## Experiment directories (`expr/`)

`expr/` is gitignored and holds one self-contained snapshot per run, named by timestamp
(`DD_MM-HH_MM`). **`./CreateExperiment <mesh.geo> <phys_file> <interface_dat_file>`** is the
normal way to start a new run: it copies `src/{cuda_cfd.cuf,interface.f90,make_vtk.f90}`, the
root `Makefile`, `tupi.slurm`, and `src/build_vtks.sh` into a fresh `expr/<timestamp>/`, renames
the given mesh/phys/interface files to `tcc_lt.geo`/`tcc_lt.phys`/`INTERFACE.dat`, and runs
`make -C expr/<timestamp> mesh interface`. Each experiment directory is then edited and built
independently (its own `PARAMETER` block, its own mesh) and is the actual thing rsynced to and
run on the cluster — treat `src/` and the root `Makefile` as the templates experiments are
stamped from, not as something built/run directly in place (except via the Makefile targets below).

## Common commands (from repo root)

- `make mesh` — `gmsh tcc_lt.geo -3 -format msh2 -o tcc_lt.msh` (expects `tcc_lt.geo` in cwd;
  experiment dirs have this, the repo root does not — run mesh-building from within an `expr/*`
  dir, or point Gmsh at a file under `mshs/` directly).
- `make interface` — compiles `src/interface.f90` with `gfortran` → `./interface`.
- `make cfd` — compiles `src/cuda_cfd.cuf` with `nvfortran` → `./cfd` (needs NVIDIA HPC SDK;
  only works where `nvfortran`/a CUDA GPU are available, e.g. the `tupi`/`pcad` cluster).
- `make run_sanitizer` — builds `cfd` then runs `compute-sanitizer` (memcheck, full leak check)
  against it, logging to `./sw/compute_sanitizer.txt`.
- `make vtk` — compiles `src/make_vtk.f90` with `gfortran` → `./vtks`.
- `make ssh` — rsyncs the whole current directory to `gppd-hpc.inf.ufrgs.br`, then SSHes into
  `pcad.inf.ufrgs.br` to `sbatch tupi.slurm`, then rsyncs matching result files
  (`tcc_lt[0-9][0-9][1-9]*`) back into `./resultados`. This touches remote shared
  infrastructure — confirm with the user before running it.
- `make tar` — tars the current directory (excluding the tarball itself) for archiving/transfer.

## Conventions worth knowing

- **`pname`**: a fixed 6-character project/case name used as a filename prefix everywhere
  (`<pname>000.cor`, `<pname>.con`, `<pname>CC.bv`, `<pname>CS01.nnn`, ...). Code validates the
  length strictly; keep it 6 characters when creating new cases.
- Comments, prompts, and variable docs in the Fortran/Python sources are in Portuguese; this is
  consistent throughout and not something to "fix."
- Field precision: the preprocessor writes `float64`-precision text (`%.16e`/`%.17e`), but the
  CUDA solver's main arrays are `real(4)` — the solver reads/operates in single precision even
  though inputs are written with double-precision formatting.
- Boundary conditions are resolved per-node by merging all physical groups a node belongs to,
  keeping the most restrictive (smallest-magnitude, or any non-free) condition per velocity
  component — see `merge_condition`/`merge_nodal_bc` in `TG2_preprocessor.py`.
- Immersed-object face normals are oriented to point into the fluid domain using the object's
  bounding-box center as a reference; this assumes the immersed body is star-convex about that
  center (true for boxes/cylinders/spheres/prisms, not guaranteed for complex shapes).
