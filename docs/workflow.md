# A snapshot-to-products workflow

From the repository root, run the example with a fresh output directory.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --threads=4 --project=. examples/mhd_workflow.jl --demo /tmp/aegis-mhd-demo
```

The public fixture has two shuffled, unequal-volume MeshBlocks. Its fluid velocity is $(3,4,0)$, its field is $(0,0,2)$, and density varies with position. The run covers 189 code volume units and contains 1,008 code mass units. The kinetic energy integral is $1008(3^2+4^2)/2=12600$ and the magnetic energy integral is $189(2^2)/2=378$. The program verifies both. Its slab probe and slice use a user field called `tracer`. The demo creates its source files in a temporary directory and recreates them on each invocation.

| Output | Contents |
| --- | --- |
| `summary.toml` | Source path, stream, stored time and cycle when available, package version, model, field inventory, input blocks and history schemas when present, selection, one native cell probe, volume, integrated energies and execution settings. |
| `mass_by_density.h5` | Gas mass per density bin and normalized probability density per dex of code density. The default physical bin edges range from $10^{-4}$ to $10^4$. Included and excluded masses show how much the chosen edges cover. |
| `slice.h5` | Native cell layer, physical axis/coordinate, transverse pixel faces and centers, field metadata and source selection. The default cut is perpendicular to $x_2$ at the selected center closest to the clipped interval midpoint. An exact distance tie selects the lower center. |

Replace `--demo` with a run directory for your own data. The command selects the last unambiguous time in a single discovered mesh stream. For a real run the default slice shows magnetic energy density. A third argument replaces it with an available stored field. The function accepts an explicit stream, region, physical bin edges, slice axis and coordinate. To slice a derived field, pass a function such as `slice_field=s -> quantity(s, NewtonianMHD(), :speed)`.

```julia
using Aegis
include("examples/mhd_workflow.jl")
using .MHDWorkflowExample: run_mhd_workflow

result = run_mhd_workflow(run_directory, new_output_directory;
    stream="run.out2", slice_field="rho", slice_axis=1, at=0.0,
    region=Region(x2=(-1,1)), density_edges=[0.01, 0.1, 1.0, 10.0, 100.0],
    workers=4, chunk_blocks=2)
result["kinetic_energy"], result["magnetic_energy"]
```

The workflow declares Newtonian fluid primitives in code units and reads `rho`, `vel1`–`vel3`, and `Bcc1`–`Bcc3` from the selected snapshot. It uses the validated Cartesian geometry for the integral, PDF and slice. A mesh-only directory works with the same model declaration. Histories and input files, when available, supply independent context. A snapshot time collision requires the caller to choose a state outside this last-snapshot example. [Quantities](quantities.md) · [Distributions and spatial products](products.md)

## Validation on the supplied runs

```sh
julia --threads=4 --project=. test/workflow_integration.jl
```

The integration suite requires all three runs. Each comparison below uses an exact unique printed history time equal to the stored final snapshot time, complete root-domain coverage, and the verified Newtonian simulation frame. In the compatible solvers, the history writer accumulates three components $\sum_i \Delta V_i m_{a,i}^2/(2\rho_i)$ and three components $\sum_i \Delta V_i B_{\mathrm{cc},a,i}^2/2$. Here $m_a=\rho v_a$ and the snapshot contains primitive velocity and cell-centered magnetic components. The sum of the three printed directional integrals provides an independent check on each energy. Tolerances allow half of the last printed digit in each component and rounding of native Float32 snapshot values. [MHD history writer](https://github.com/dfielding14/athena-1/blob/d2aa0b0ffc3b6b215721693cf82946cfa5582616/src/outputs/history.cpp) · [MHD-PIC history writer](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/history.cpp)

| Final mesh | Code time | Kinetic energy, snapshot / history | Magnetic energy, snapshot / history | Slice |
| --- | ---: | ---: | ---: | --- |
| 128³ MHD | 10 | 921.872316 / 921.872 | 121.688813 / 121.6888 | `ohmic_heating` |
| 256³ MHD | 10 | 863.196751 / 863.198 | 104.253164 / 104.2531 | `ohmic_heating` |
| 1024 × 4096 × 1 MHD-PIC | 600 | 2274.630795 / 2274.626 | 304.846329 / 304.846 | magnetic energy density |

All six energy comparisons passed. Every density histogram captured the full gas mass to the solver history's printed precision. The two MHD slices use the original problem-generator field discovered in the snapshot. Particle energies and custom history sums have separate definitions.

## Timing and memory

```sh
run=test/test_data/3D_Athena/256_eta5e-3
julia -t 4 --project=. benchmark/workflow.jl "$run" 1 1
julia -t 4 --project=. benchmark/workflow.jl "$run" 4 1
julia -t 4 --project=. benchmark/workflow.jl "$run" 4 4
```

Each configuration starts a fresh Julia process, warms the code on the synthetic example, then runs the real workflow four times with identical selections and output contents. Products go to new temporary directories. The timed call includes stream discovery, input/history parsing, geometry, field reads, both energy reductions, the mass-weighted PDF, the magnetic-energy slice and all writes. Garbage collection runs before each trial. Compiling, decompression and filesystem caching can affect first-access timing. Cached writes are measured through file close and rename. End-to-end throughput divides the number of distinct selected cells by total elapsed time. Each cell participates in multiple calculations.

A single task owns the HDF5 file. The native reductions can send chunks of materialized block data to compute tasks in the same process. At most $W K$ derived blocks are in flight in a wave with $W$ workers and $K$ blocks per task. A $32^3$ Float64 energy-density block occupies 256 KiB; construction of a kinetic block also reads four 128 KiB Float32 primitive arrays. The distribution visits blocks serially. Additional memory includes the Cartesian geometry, field metadata, parsed histories, histogram bins, the transverse slice, storage-library buffers and the Julia runtime. High-water RSS includes the entire process and all threads. [Native benchmark details](performance.md)

**Measured 256³ example.** A 32 GB Apple M2 Max Mac Studio (12 CPU cores) ran Julia 1.12.0 and HDF5.jl 0.17.3 with four Julia compute threads in all configurations. The file contains 512 blocks of 32³ cells, with 16,777,216 distinct cells in the selection. Times and throughput report the median of three warm trials. RSS is the greatest process high-water reading through the last trial. Allocations are cumulative Julia allocations during one warm trial.

| Workers | Blocks per task | Synthetic warm-up (s) | First real (s) | Warm end-to-end (s) | Million cells/s | Allocated (MiB) | Peak RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 3.40 | 1.39 | 1.202 | 13.96 | 787.0 | 628.5 |
| 4 | 1 | 3.22 | 1.13 | 0.966 | 17.36 | 789.5 | 648.6 |
| 4 | 4 | 3.19 | 1.18 | 0.952 | 17.63 | 788.4 | 655.8 |

Across these settings, the integrated kinetic and magnetic energies agreed to the reported digits. Only the reductions parallelize their arithmetic. Discovery, reads, the PDF, the slice and export contribute to all three wall times. These measurements characterize the cached local file and the stated worker settings.
