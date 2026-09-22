# HYDRUS-1D in Rust

A modern port of the HYDRUS-1D computational engine (Šimůnek, van Genuchten & Šejna) to Rust,
with a new cross-platform GUI (egui) that builds on **Linux and Windows**.

```
crates/
  hydrus-core   numerical engine (water flow, root uptake, solute, heat) – no I/O, no GUI
  hydrus-io     legacy HYDRUS-1D project reader (Selector.in / Profile.dat / Atmosph.in),
                native .h1dr (JSON) projects, CSV export
  hydrus-cli    hydrus1d-cli – run a project headless
  hydrus-gui    hydrus1d – the desktop application
validation/     scripts used to compare against the original Fortran code
```

## Build

```
cargo build --release            # needs Rust >= 1.80
./target/release/hydrus1d        # GUI
./target/release/hydrus1d-cli <legacy-project-dir | project.h1dr> -o out/
```
Linux: `libxkbcommon`, `libwayland`, `libx11`, OpenGL are needed at run time (already present on desktop distributions).
Windows: nothing extra. `.github/workflows/build.yml` builds and packages both platforms (and attaches
them to a release when you push a `v*` tag).

## What is ported (from the supplied Fortran sources)

| Area | Status |
|---|---|
| Richards equation, Picard iteration, time-step control, all BC switching logic (atmospheric with/without surface layer, constant/variable head & flux, free drainage, seepage face, GWL-flux, 5 drain types) | ✔ ported and compared to Fortran |
| Hydraulic models: van Genuchten, modified VG, Brooks–Corey, VG with air entry, Kosugi, Durner; property look-up tables | ✔ |
Hysteresis: Kool–Parker style scanning curves (iHyst 1 and 2); Lenhard air-entrapment hysteresis (iHyst 3) | ✔ (Kool–Parker option 1 compared to Fortran; option 2 and Lenhard iHyst 3 ported but not yet cross-validated) |
| Root water uptake (Feddes, S-shaped, compensated uptake, root growth, solute stress) | ✔ |
| Solute transport: Galerkin / upstream weighting, non-linear sorption, decay chains, gas phase, two-site and mobile-immobile non-equilibrium, mass-balance | ✔ compared to Fortran |
| Heat transport (conduction, convection, sinusoidal / atmospheric top temperature) | ✔ compared to Fortran |
| Legacy import (Selector.in v3/v4, Profile.dat, Atmosph.in) | ✔ |
| **Not ported yet** | vapor flow, Meteo/Penman–Monteith BC, snow, interception/LAI, dual-porosity & dual-permeability, virus/colloid/filtration, temperature- and moisture-dependent reaction rates, active root solute uptake, inverse module, legacy *writer*, `Hysteresis.in`/`Options.in`/`MoistDep.in` reading |

Unsupported options in an imported project produce a clear error instead of silently different results.

## Validation

`validation/build_reference.sh` compiles the original Fortran with gfortran (patching only the MS-Fortran
system calls), the `gen*.py` scripts create identical input directories, and `cmp*.py` compare results.
On the test cases used during the port (infiltration/evaporation with layered soil and root stress, hysteresis,
constant head, seepage face, GWL flux, drains; two-solute chain with non-linear sorption and two-site kinetics;
heat transport) cumulative fluxes and profiles agree to 4–5 significant digits until the two codes take
(chaotically) different adaptive time steps. The original uses single precision, the port uses f64, so
bit-identical results are not expected. Water-balance errors of the test problem are identical to the Fortran ones.

Run the tests with `cargo test --release -p hydrus-core -p hydrus-io`.

## GUI status

Tested on Linux (Xvfb + software OpenGL): project import, all main pages, running, profile / flux plots.
The Windows build is produced by CI and has not been run by the author. `hydrus1d <folder|file.h1dr>` opens a
project directly from the command line.

## Notes on conventions

* Units are user units, as in HYDRUS (no internal conversion), except where the original converts (`xConv`, `tConv`).
* Projects list nodes top → bottom (depth positive downward); the solver stores them bottom → top like the Fortran.
* Constant-head boundaries take the value from the initial condition at the boundary node (as HYDRUS does).
