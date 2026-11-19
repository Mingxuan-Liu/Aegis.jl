# Aegis.jl — Initial MHD implementation roadmap

## Scope and first useful result

Begin with the MHD analysis needed across multiple projects. Particle analysis, particle-output conventions, and ML training are deferred. This is a current priority decision, not a permanent package restriction; the lasting development principles are in `AGENTS.md`.

The first milestone is a complete, documented workflow that reads a run's history and selected mesh data, interprets the fields correctly, computes volume integrals and weighted statistics, produces a basic slice or projection, and saves results with provenance. It must work within a stated memory budget and have an independently validated serial path plus one measured parallel path.

**Current status.** Items 1–4 are implemented for the supported native output layouts and complete, unrefined Cartesian reductions. The public example and tests run on generated data, and the integration suite covers the two MHD runs and the fluid output of the MHD-PIC run. Item 5, fundamental derived MHD quantities, is next. The slice, projection and export portions of the first milestone are scheduled under item 6. [User guide](README.md) · [Input contract](docs/input-contract.md) · [Benchmark](docs/performance.md)

## Priority 1 — Establish the data contract and a working foundation

### 1. Inspect representative outputs and declare supported cases

- Obtain small representative `.athdf`, `.hst`, and input-parameter files from the actual Athena++ branch used for these projects. Record the solver version or commit when available.
- Verify output names, dataset axes, primitive and conserved fields, magnetic normalization, units, equation of state, time information, mesh geometry, refinement metadata, and field centering against the writer and a trusted reader.
- Define explicit capability checks. A reasonable first supported case is nonrelativistic MHD on Cartesian meshes; confirm this against actual inputs. Support a narrow case correctly and reject unsupported cases clearly.
- Include non-cubic, multi-block examples and singleton dimensions. Preserve refinement information in the data model even if some initial diagnostics support only uniform grids.

**Acceptance:** a short input specification and small, provenance-labelled fixtures make physical and storage conventions unambiguous. No routine needs to guess a unit conversion, boundary condition, axis, or required EOS parameter.

### 2. Create the Julia package and its validation scaffold

- Establish a normal Julia package environment, public module, dependency compatibility, automated tests, and documentation structure.
- Give readers, physical data representation, scientific operations, and user workflows distinct responsibilities. Keep plotting and heavyweight optional backends from becoming unnecessary core requirements.
- Establish small deterministic fixtures, a way to run independent reference checks, and a reproducible benchmark environment.
- Add a user-facing README as the first workflow becomes executable; only document installation and API examples that have actually been tested.

**Acceptance:** a clean checkout can instantiate the environment and run its existing tests. The initial interface has explicit inputs and does not depend on private paths or an interactive session.

### 3. Implement metadata-aware readers

- Read `.athdf` metadata and requested fields without requiring full-domain reconstruction. Support native mesh-block iteration and bounded-size reads; distinguish physical dimensions from storage-array order.
- Read `.hst` headers and named columns, including user-defined quantities. Make duplicate times, restart overlaps, repeated headers, and out-of-order records visible, with a documented policy and access to raw records.
- Discover and order a run's outputs using validated time metadata. Avoid treating filenames alone as evidence of temporal ordering or synchronization.
- Read the needed simulation parameters or accept an explicit configuration when parameters are not embedded in outputs.

**Acceptance:** selected values, coordinates, shapes, and times agree with independently constructed fixtures and Athena++'s reference reader. Partial reads have a measured, bounded memory cost consistent with the requested data, storage chunks, and necessary metadata, without materializing the full snapshot.

## Priority 2 — Deliver reusable MHD calculations

### 4. Implement geometry-aware selection and reductions

- Represent active cells, cell volumes, spatial bounds, and field centering. Separate selection from calculation.
- Provide integrals, extrema, moments, and explicitly volume- or mass-weighted means. Make empty selections and invalid-value handling explicit.
- Implement chunked reductions whose partial results can be combined correctly. A parallel mean must combine the weighted numerator and denominator, not average worker means blindly.
- For supported refined meshes, reduce over valid physical coverage once. Specify the behavior of selections that intersect cells and the distinction between native-grid and resampled quantities.

**Acceptance:** constant and manufactured fields reproduce known integrals on supported geometries, including unequal cell sizes where claimed. Chunk size, block order, and worker partitioning do not change the scientific result beyond the justified tolerance.

### 5. Add fundamental derived MHD quantities

- Provide the quantities justified by available fields and the chosen model: velocity and magnetic magnitudes, kinetic and magnetic energy densities, magnetic pressure, plasma beta, Alfvén speed, and relevant characteristic speeds.
- Add thermal energy, sound speed, temperature, and Mach numbers only with explicit EOS and conversion information. Missing pressure, composition, or thermodynamic parameters must not be silently guessed.
- Keep field dependencies visible so a calculation reads only what it needs. Distinguish local densities from volume-integrated totals.
- Define behavior for zero magnetic fields, invalid density or pressure, non-finite inputs, and quantities outside their physical domain. Do not hide problems with unexplained clipping or floors.

**Acceptance:** independent uniform-state examples and limiting cases verify each formula, its units, and its magnetic/EOS convention. Tests cover both valid results and defined failure behavior.

### 6. Add distributions, conditional statistics, and basic spatial products

- Implement weighted histograms/PDFs and conditional reductions using the same selection and reduction machinery. Specify whether an output is a bin-integrated amount or a density per unit linear/logarithmic coordinate.
- Provide reusable joint distributions and phase-conditioned statistics without embedding a particular temperature window or research interpretation.
- Produce slices and line-of-sight integrals or weighted projections with explicit axes, coordinates, path lengths, and units. Keep the numerical result separate from plotting style.
- Export compact result arrays/tables with labels, units, and provenance; use an explicit, documented output convention.

**Acceptance:** distribution normalization closes under the stated measure; analytical slab examples verify projections; and users can change selections, weighting, and fields without rewriting a calculation.

### 7. Complete one user workflow and measure its scaling

- Demonstrate the path from a run directory to history inspection, selected-field access, energy totals, one weighted distribution, and one spatial product.
- Add a safe parallel path over independent blocks, chunks, or snapshots. Choose the decomposition from measured work and I/O behavior; control memory ownership and concurrent reads.
- Compare serial and parallel results and measure representative runtime, allocations or peak memory, throughput, and scaling. Separate compilation and storage effects from kernel timing.
- Have independent reviewers challenge the physical definitions, mesh treatment, tests, and performance claims. Incorporate the resulting durable fixes and regression cases.
- Write the README around this executable workflow and the relationship between its components. Remove temporary working artifacts while keeping validation and benchmark assets.

**Acceptance:** a researcher can reproduce the example from documented inputs. At least one representative workload has a recorded correctness check and measured serial/parallel performance, with limitations stated honestly.

## Priority 3 — Extend the verified MHD foundation

### 8. Add differential operators

- Add gradients, divergence, curl, vorticity, and magnetic gradients for explicitly supported meshes and boundary conditions.
- Handle centering, block interfaces, boundary stencils, and refinement transitions deliberately. Do not compute derivatives independently inside blocks while ignoring neighboring data.
- Distinguish derivatives reconstructed from output fields from native discrete solver operators. In particular, divergence from cell-centered magnetic fields is not automatically the solver's face-based divergence diagnostic.
- Add current or dissipative-heating estimates only after confirming coefficient definitions, magnetic conventions, constitutive assumptions, and whether the output supports the required reconstruction.

**Acceptance:** constant, linear, and smooth periodic manufactured fields verify identities and convergence. Block boundaries do not introduce artificial discontinuities in supported configurations; unsupported boundaries or refinement cases fail clearly.

### 9. Add Fourier and spectral diagnostics

- Start with uniform Cartesian periodic data, or an explicitly documented resampling step. Do not label a transform of concatenated mesh blocks as a physical-domain spectrum.
- Define the Fourier normalization, physical wavenumbers, mean/DC treatment, shell/bin measures, and whether output is integrated power per bin or spectral density.
- Account correctly for a real-valued transform's conjugate modes and exceptional zero/Nyquist modes. Verify behavior on non-cubic grids and unequal domain lengths.
- Validate vector-field and kinetic/magnetic spectra against the exact real-space quantity each definition is intended to represent. Document the choice of density weighting for kinetic spectra.
- Estimate peak FFT memory and parallel-library settings. Add distributed transforms only when representative workloads demonstrate the need.

**Acceptance:** a single Fourier mode appears at the expected physical wavenumber, and integrated spectral power agrees with the corresponding real-space norm under the documented normalization.

## How testing fits the schedule

Scientific validation, adversarial review, documentation, and memory/performance checks accompany each relevant item. They are not a final phase after feature development. Keep rapid deterministic tests in routine CI and larger integration/scaling checks in a documented reproducible suite.

## Suggested repository responsibilities

| Location | Purpose and relationship |
| --- | --- |
| `src/` | Reusable readers, physical data interfaces, and numerical operations; imported by examples and downstream research code. |
| `test/` | Independent fixtures, analytic checks, regressions, and integration tests for supported behavior. |
| `examples/` | Small runnable user workflows that compose the public package interface. |
| `docs/` | Scientific definitions, conventions, assumptions, and usage explanations. |
| `benchmark/` | Repeatable workloads for measuring time, memory, I/O, and scaling. |
| `README.md` | The user entry point: purpose, installation, a working example, and guidance on which components to use. |

Create directories when they have useful content. Project-specific publication scripts belong in downstream projects or clearly scoped examples, with shared calculations imported from the package.

## Reference starting points

- [Athena++ reference readers](https://github.com/PrincetonUniversity/athena/blob/main/vis/python/athena_read.py): useful for input interpretation and independent comparisons; confirm against the actual solver branch.
- [Julia performance guidance](https://docs.julialang.org/en/v1/manual/performance-tips/): computational structure, allocations, memory access, and benchmarking.
- [Julia threading guidance](https://docs.julialang.org/en/v1/manual/multi-threading/): task ownership and race-free parallel work.
- [HDF5.jl parallel I/O](https://juliaio.github.io/HDF5.jl/stable/mpi/): backend requirements for MPI-based I/O; distinct from simply using Julia threads.

Specific scientific diagnostics require their own literature-backed definitions under `AGENTS.md`.
