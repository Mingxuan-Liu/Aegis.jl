using Aegis

# Write a small run with known answers into a temporary directory.
include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

mktempdir() do dir
    files = example_run(dir)
    input = read_athinput(files.input)
    @show input_blocks(input)
    @show parameter(input, "problem", "tracer_scale", Float64)

    history = read_history(files.history)
    @show history_column(history, "time")
    @show keys(discover_snapshots(dir))

    open_snapshot(files.snapshot) do s
        @show s.time s.blocksize [(f.dataset, f.name) for f in s.fields]
        tracer = interpret_field(s, "tracer"; alias="sample", units="code tracer unit")
        @show read_field(s, tracer, 1; ranges=(1:2, 1:1, 1:1))
        region = Region(x1=(1, 6))
        mesh = CartesianMesh(s)
        result = summarize(mesh, tracer; region)
        by_mass = summarize(mesh, tracer; region, density="rho")
        @show result.volume result.integral result.mean
        @assert result.volume == 45 && result.integral == 117
        @assert result.mean ≈ 2.6
        @show by_mass.mean
        @assert by_mass.mean ≈ 8/3
    end
end
