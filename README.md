# Aegis.jl

Aegis is a Julia package for analyzing the fluid in Athena++ MHD and MHD-PIC simulations. Use it to compute Newtonian energies, weighted distributions, slices and projections on native Cartesian cells. It can also reconstruct gradients across MeshBlocks and measure three-dimensional Fourier spectra on uniform periodic grids.

## Analyze a run

Use Julia 1.10 or newer. From the repository root, run the example with a new output directory.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --threads=4 --project=. examples/mhd_workflow.jl --demo /tmp/aegis-mhd-demo
```

The example analyzes an uneven two-block grid in shuffled storage order. From the gas density, velocity and cell-centered magnetic field it recovers total kinetic energy **12,600**, magnetic energy **378**, and mass **1,008**, all in code units. It writes three files.

| File | What to inspect |
| --- | --- |
| `summary.toml` | Snapshot identity, mesh and field catalog, selection, integrated energies, and any accompanying input and history files |
| `mass_by_density.h5` | Mass in each density bin, the normalized PDF per dex, and mass outside the bin edges |
| `slice.h5` | A cell-centered slice with its transverse coordinates, units when declared, and selection metadata |

To analyze your own Newtonian primitive output, replace `--demo` with the run directory. The optional final argument selects a stored field for the slice. The following example uses a column produced by the filament-heating problem generator.

```sh
julia --threads=4 --project=. examples/mhd_workflow.jl test/test_data/3D_Athena/128_eta5e-3 /tmp/aegis-filament-demo ohmic_heating
```

The workflow selects the latest unambiguous state of a single stream and can operate on mesh files alone. Its Julia entry point accepts a stream, region, density-bin edges and slice position. [Workflow and measured scaling](docs/workflow.md)

## Choose an analysis

| Scientific question | Julia interface and assumptions |
| --- | --- |
| How much magnetic or kinetic energy is in a region? | Pass a `quantity` from a declared `NewtonianMHD` model to `summarize`. Thermodynamic quantities require an appropriate EOS and physical scales where applicable. [Quantities](docs/quantities.md) |
| How does mass or volume depend on one or two fields? | `distribution` and `conditional_stats` accept native, derived, or caller-defined fields. `slice_plane` and `project` use physical cell lengths. [Distributions and spatial products](docs/products.md) |
| What is the vorticity across a MeshBlock boundary? | Build a `DifferentialGrid` with explicit domain boundary rules, then read the components returned by `curl`. [Differential operators](docs/differentials.md) |
| How is magnetic energy distributed over spatial scales? | `shell_spectrum` takes three cell-centered magnetic fields and confirmed periodic boundaries. It reports energy per shell and per unit angular wavenumber, including a real-space normalization check. [Spectra](docs/spectra.md) · [Final-snapshot figure comparison](validation/spectra/README.md) |

## Inspect individual inputs and cells

Run `julia --project=. examples/native_cells.jl` for a smaller example of the reader and reduction interfaces. The example creates a temporary two-block run, reads its custom input parameter and history, inspects the snapshot fields, and measures a passive scalar called `tracer`. The blocks are stored out of spatial order. Its six cells along $x_1$ have widths 1 through 6 and tracer values 1 through 6. The selection $1 \leq x_1 < 6$ includes two cell centers. Their widths are 2 and 3, and the transverse area is 9. The selected volume is 45, the tracer integral is $9(2\times2+3\times3)=117$, and the volume-weighted mean is 2.6.

For your own run, start by inspecting the output streams and the input file. Each history row retains its original header and line number.

```julia
using Aegis

streams = discover_snapshots(run_directory)
keys(streams)                                      # one series per output stream
input = read_athinput(input_path)
input_blocks(input)
parameter_names(input, "problem")
parameter(input, "problem", "tracer_scale", Float64)  # for the example input

h = read_history(history_path)
h.schemas[1].names
history_column(h, "time")                          # original rows, in file order
history_column(h, "time"; rows=restart_rows(h))    # explicit restart policy
```

A snapshot is self-describing. You can open one without a history, input deck, pressure field or particle file. A field descriptor identifies one variable slot in one dataset and snapshot. The native field reader returns an ordinary Julia array in $(x_1,x_2,x_3)$ order and preserves the stored precision. A single-cell range keeps its axis.

```julia
open_snapshot(snapshot_path) do s
    [(f.dataset, f.index, f.name, f.centering, f.precision) for f in s.fields]
    s.time, s.cycle, s.coordinates, s.rootsize, s.blocksize
    first_block = block_coordinates(s, 1)
    slab = read_field(s, "tracer", 1; ranges=(1:2, 1:1, 1:1))

    # Keep lab-specific meanings and units alongside their file identities.
    scalar = interpret_field(s, "tracer"; alias="sample", units="code tracer unit")
    mesh = CartesianMesh(s)
    region = Region(x1=(1, 6))
    by_volume = summarize(mesh, scalar; region)
    by_mass = summarize(mesh, scalar; region, density="rho")
    (; by_volume.integral, volume_mean=by_volume.mean, mass_mean=by_mass.mean)
end
```

Substitute fields and bounds from your output. Use a unique name or a pair such as `("hydro", 2)` giving the dataset and one-based variable position. `interpret_field` records what you supply. Its labels do not convert the stored numbers. For mass weighting, pass a field you have verified represents gas mass density in the chosen units and frame.

## What a reduction measures

`CartesianMesh(s)` verifies level-zero block coverage and the stored coordinates. Reuse it for multiple fields in the same open snapshot. `selectedblocks(mesh; region)` exposes the block indices and ranges for custom analyses. `summarize` accepts either the open snapshot or its verified mesh. Set `workers` and `chunk_blocks` to control threaded arithmetic and memory use.

| Result | Meaning |
| --- | --- |
| `count`, `volume` | Number and total coordinate volume of selected cells |
| `integral` | Sum of the stored field times cell volume, including when a density is supplied |
| `weight` | Total volume or total mass used in the moments |
| `mean`, `variance` | Weighted mean and weighted population variance |
| `minimum`, `maximum` | Extremal selected field values |
| `provenance` | Source, selection, field and weighting information, and execution settings |

Region bounds are lower-inclusive and upper-exclusive in each coordinate. Each selected center contributes its whole cell. Unequal face spacing and singleton dimensions use their recorded lengths. An integral has units of field times coordinate volume, the mean has field units, and the variance has squared field units. The supported reduction case is cell-centered output on a complete unrefined Cartesian mesh of active cells. The raw reader preserves native block levels and logical locations for inspection of other compatible layouts.

## Scientific conventions and validation

The [input contract](docs/input-contract.md) describes the inspected writers, the two observed layouts, physical conventions, support limits and the three supplied runs. The [performance notes](docs/performance.md) give benchmark commands and memory measurements. See [the roadmap](MHD_ROADMAP.md) for later priorities.

Run the public tests from any checkout. They generate their own inputs locally.

```sh
julia --threads=4 --project=. -e 'using Pkg; Pkg.test()'
```

The separate integration test reads the supplied run directories under `test/test_data/`. It requires Python with NumPy and h5py for independent HDF5 value and integral checks. Those dependencies are used only for validation.

```sh
julia --threads=4 --project=. test/integration.jl
julia --threads=4 --project=. test/workflow_integration.jl
julia --threads=4 --project=. test/differentials_integration.jl
```

To also compare selected cells with Athena++ itself, set `AEGIS_ATHENA_READER` to a compatible local `vis/python/athena_read.py` before running `test/integration.jl`. `PYTHON` can select the interpreter. A positional argument overrides the private data root. The suite requires all three run directories.

For a scientific citation, use Stone et al. (2020) for the Athena++ MHD framework and Sun & Bai (2023) for the MHD-PIC method. Record the Aegis revision, solver provenance, source snapshot and any user-supplied units alongside analysis results. [Stone et al.](https://doi.org/10.3847/1538-4365/ab929b) · [Sun & Bai](https://arxiv.org/abs/2304.10568)
