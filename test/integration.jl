# Optional, read-only checks for the three supplied runs. Requires h5py + NumPy.
# julia -t 4 --project=. test/integration.jl [parent-of-the-three-run-directories]
using Aegis
using Test
using TOML

const REAL_RUNS = [
    (directory="3D_Athena/128_eta5e-3", grid=(128, 128, 128), quantity="ohmic_heating",
     parameter=0.005, custom_history="Mean Ms", column=14, lasttime=10.0),
    (directory="3D_Athena/256_eta5e-3", grid=(256, 256, 256), quantity="ohmic_heating",
     parameter=0.005, custom_history="Mean Ms", column=14, lasttime=10.0),
    (directory="2D_MHD_PIC/run424", grid=(1024, 4096, 1), quantity="vel1",
     parameter=0.0, custom_history="Pideal_y", column=29, lasttime=600.0),
]

function check_run(parent, case)
    dir = joinpath(parent, case.directory)
    streams = discover_snapshots(dir)
    isempty(streams) && error("No .athdf snapshots found in $dir")
    # All three supplied runs have one mesh output stream. A newly added
    # unsynchronized stream needs its own explicit selection.
    records = only(values(streams))
    @test issorted([r.time for r in records])
    record = last(records)
    @test record.time == case.lasttime
    files = readdir(dir; join=true)
    inputs = filter(f -> startswith(basename(f), "athinput"), files)
    histories = filter(f -> endswith(f, ".hst"), files)
    input = read_athinput(only(inputs))
    @test Tuple(parameter(input, "mesh", "nx$d", Int) for d in 1:3) == case.grid
    @test parameter(input, "problem", "eta_ohm", Float64) == case.parameter
    history = read_history(only(histories))
    custom_values = history_column(history, case.custom_history)
    @test length(custom_values) == length(history.rows)
    last_raw_row = open(only(histories)) do io
        first(Iterators.reverse(eachline(io)))
    end
    @test last(custom_values) == parse(Float64, split(last_raw_row)[case.column])
    history_times = history_column(history, "time")
    history_index = only(findall(==(record.time), history_times))

    python = get(ENV, "PYTHON", "python3")
    script = joinpath(@__DIR__, "reference.py")
    cmd = `$python -B $script $(record.path) $(case.quantity)`
    reader = get(ENV, "AEGIS_ATHENA_READER", "")
    isempty(reader) || (cmd = `$cmd --athena-reader $reader`)
    ref = TOML.parse(read(cmd, String))

    open_snapshot(record.path) do s
        @test s.rootsize == case.grid
        @test collect(s.blocksize) == ref["blocksize"]
        @test (s.time, s.cycle) == (ref["time"], ref["cycle"])
        @test [f.name for f in s.fields] == ref["names"]
        @test [f.dataset for f in s.fields] == ref["datasets"]
        block = ref["block"]
        coords = block_coordinates(s, block)
        inds = ref["indices"]
        for d in 1:3
            i = inds[d]
            @test coords.faces[d][i:i+1] == ref["faces"][d]
            @test coords.centers[d][i] == ref["centers"][d]
        end
        ranges = ntuple(d -> inds[d]:inds[d], 3)
        for (f, expected) in zip(s.fields, ref["samples"])
            slab = read_field(s, (f.dataset, f.index), block; ranges)
            @test size(slab) == (1, 1, 1)
            @test eltype(slab) == f.precision == Float32
            @test only(slab) == expected
        end

        mesh = CartesianMesh(s)
        x, y, z = Tuple.(ref["roi"])
        local_q = summarize(mesh, case.quantity; region=Region(x1=x, x2=y, x3=z))
        @test local_q.count == ref["roi_count"]
        @test isapprox(local_q.volume, ref["roi_volume"]; rtol=1e-14)
        @test isapprox(local_q.integral, ref["roi_integral"]; rtol=1e-14)
        rho = summarize(mesh, "rho")
        q = summarize(mesh, case.quantity; density="rho", workers=min(4, Threads.nthreads()),
                      chunk_blocks=2)
        @test q.count == rho.count == prod(case.grid)
        # Accumulators use Float64; per-block partials cover up to 32^3 cells.
        # Allow roundoff in different summation orders, scaled by field range
        # for means (the signed PIC velocity mean is close to zero).
        @test isapprox(rho.volume, ref["volume"]; rtol=2e-12)
        @test isapprox(rho.integral, ref["mass"]; rtol=2e-12)
        @test isapprox(q.weight, ref["mass"]; rtol=2e-12)
        @test q.minimum == ref["minimum"]
        @test q.maximum == ref["maximum"]
        scale = max(abs(q.minimum), abs(q.maximum))
        @test isapprox(q.mean, ref["mass_mean"]; rtol=2e-12, atol=2e-13*scale)

        # Solver's standard 'mass' at this exact snapshot time covers the same
        # full domain. Histories here print six significant figures. Allow half
        # the last printed digit and the separate rounding of archived Float32 data.
        solver_mass = history_column(history, "mass")[history_index]
        half_digit = 0.5 * 10.0^(floor(Int, log10(abs(solver_mass))) - 5)
        @test abs(rho.integral - solver_mass) <= half_digit + 3eps(Float32)*abs(solver_mass)
        println(case.directory, "  t=", s.time, "  blocks=", s.nblocks, "  fields=",
                length(s.fields), "  history rows=", length(history.rows),
                "  volume=", rho.volume, "  mass=", rho.integral,
                "  mass-weighted ", case.quantity, "=", q.mean)
    end
end

length(ARGS) <= 1 || error("usage: test/integration.jl [test_data_directory]")
root = isempty(ARGS) ? joinpath(@__DIR__, "test_data") : only(ARGS)
for case in REAL_RUNS
    isdir(joinpath(root, case.directory)) || error(
        "Expected private run $(case.directory) is missing under $root. Restore it or supply its data root.")
end
@testset "Supplied Athena++ runs" begin
    for case in REAL_RUNS
        @testset "$(case.directory)" begin
            check_run(root, case)
        end
    end
end
