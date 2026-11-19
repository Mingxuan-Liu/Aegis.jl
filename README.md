# Aegis.jl

Aegis reads native Athena++ mesh output for researchers working with MHD and the fluid component of MHD-PIC simulations. It can inspect a run, read a rectangular slab of a named field from one MeshBlock, and compute integrals and weighted statistics over a Cartesian region. A simulation can use its own problem-generator fields, history columns and input parameters.

## Run a first analysis

Use Julia 1.10 or newer. From the repository root, instantiate the package environment and run the example.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/native_cells.jl
```

The example creates a temporary two-block run, reads its custom input parameter and history, inspects the snapshot fields, and measures a passive scalar called `tracer`. The blocks are stored out of spatial order. Its six cells along $x_1$ have widths 1 through 6 and tracer values 1 through 6. The selection $1 \leq x_1 < 6$ includes two cell centers. Their widths are 2 and 3, and the transverse area is 9. The selected volume is 45, the tracer integral is $9(2\times2+3\times3)=117$, and the volume-weighted mean is 2.6.

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
```

To also compare selected cells with Athena++ itself, set `AEGIS_ATHENA_READER` to a compatible local `vis/python/athena_read.py` before running that command. `PYTHON` can select the interpreter. A positional argument overrides the private data root. The suite requires all three run directories.

For a scientific citation, use Stone et al. (2020) for the Athena++ MHD framework and Sun & Bai (2023) for the MHD-PIC method. Record the Aegis revision, solver provenance, source snapshot and any user-supplied units alongside analysis results. [Stone et al.](https://doi.org/10.3847/1538-4365/ab929b) · [Sun & Bai](https://arxiv.org/abs/2304.10568)
