module MHDWorkflowExample

using Aegis
using HDF5
using TOML

include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

# Nonuniform cells with unequal block volumes, v=(3,4,0), B=(0,0,2).
function write_demo(dir)
    path = joinpath(dir, "flow.out2.00001.athdf")
    fields = ("rho", "vel1", "vel2", "vel3", "Bcc1", "Bcc2", "Bcc3", "tracer")
    write_fixture(path; fields)
    h5open(path, "r+") do h
        # The fixture stores rho in "gas" and the other seven fields in "extras".
        for (slot, value) in enumerate((3, 4, 0, 0, 0, 2))
            h["extras"][:,:,:,:,slot] = fill(Float32(value), 3,2,2,2)
        end
    end
    write(joinpath(dir, "athinput.flow"), "<job>\nproblem_id = flow\n<problem>\ntracer_scale = 2.5\n")
    write(joinpath(dir, "flow.hst"),
        "# Athena++ history data\n# [1]=time [2]=mass [3]=1-KE [4]=2-KE [5]=3-KE [6]=1-ME [7]=2-ME [8]=3-ME\n1.25 1008 4536 8064 0 0 0 378\n")
    return path
end

function last_snapshot(dir, requested_stream)
    streams = discover_snapshots(dir)
    if isnothing(requested_stream)
        length(streams) == 1 || throw(ArgumentError(
            "Expected one snapshot stream. Available streams: $(collect(keys(streams))). Pass stream explicitly."))
        stream = only(keys(streams))
    else
        stream = String(requested_stream)
        haskey(streams, stream) || throw(ArgumentError("No snapshot stream named $stream"))
    end
    records = streams[stream]
    last_time = last(records).time
    count(r -> r.time == last_time, records) == 1 ||
        throw(ArgumentError("Stream $stream contains ambiguous final states at time $last_time."))
    return stream, last(records)
end

function inspect_run(dir)
    files = readdir(dir; join=true)
    histories = map(filter(p -> isfile(p) && endswith(p, ".hst"), files)) do p
        h = read_history(p)
        Dict("path" => abspath(p), "rows" => length(h.rows),
             "schema_columns" => [s.names for s in h.schemas])
    end
    inputs = map(filter(p -> isfile(p) && startswith(basename(p), "athinput"), files)) do p
        a = read_athinput(p)
        Dict("path" => abspath(p), "blocks" => input_blocks(a))
    end
    return (; histories, inputs)
end

# Select an actual native layer closest to the midpoint. A geometric midpoint
# alone can fall in a cell whose center is outside a nonuniform-grid region.
function midpoint_layer(mesh, s, region, axis)
    lo = max(region.bounds[axis][1], s.rootgrid[axis][1])
    hi = min(region.bounds[axis][2], s.rootgrid[axis][2])
    mid = lo/2 + hi/2
    best, distance = Inf, Inf
    visited = Set{Int}()
    for b in selectedblocks(mesh; region)
        segment = s.locations[axis,b.block]
        segment in visited && continue
        push!(visited,segment)
        centers = block_coordinates(s,b.block).centers[axis]
        for i in b.ranges[axis]
            x = Float64(centers[i])
            delta = abs(x-mid)
            if delta < distance || (delta == distance && x < best)
                best, distance = x, delta
            end
        end
    end
    return best
end

"""Analyze the last unambiguous state of a stream. Write results to a fresh output directory."""
function run_mhd_workflow(directory, output; stream=nothing, slice_field=nothing,
                          region=Region(), slice_axis=2, at=nothing,
                          density_edges=10 .^ range(-4.0, 4.0; length=49),
                          workers=1, chunk_blocks=1)
    dir = abspath(directory)
    dest = abspath(output)
    ispath(dest) && throw(ArgumentError("Choose a new output directory: $dest"))
    1 <= slice_axis <= 3 || throw(ArgumentError("slice_axis must be 1, 2, or 3"))
    chosen, record = last_snapshot(dir, stream)
    companions = inspect_run(dir)
    mkpath(dirname(dest))
    # Publish all three files together once every calculation and write succeeds.
    return mktempdir(dirname(dest)) do tmp
        stage = mkdir(joinpath(tmp, "products"))
        result = open_snapshot(record.path) do s
            mesh = CartesianMesh(s)
            gas = NewtonianMHD()
            kinetic = quantity(s, gas, :kinetic_energy_density)
            magnetic = quantity(s, gas, :magnetic_energy_density)
            rho = interpret_field(s, "rho"; units="code mass / code length^3")
            field = if isnothing(slice_field)
                magnetic
            elseif slice_field isa Function
                slice_field(s)
            else
                slice_field
            end
            cells = selectedblocks(mesh; region)
            isempty(cells) && throw(ArgumentError("The selection contains no native cells"))
            probe = first(cells)
            cell = ntuple(d -> first(probe.ranges[d]):first(probe.ranges[d]), 3)
            example_value = only(read_field(s, field, probe.block; ranges=cell))
            ke = summarize(mesh, kinetic; region, workers, chunk_blocks)
            me = summarize(mesh, magnetic; region, workers, chunk_blocks)
            hist = distribution(mesh, BinAxis(rho, density_edges; measure=:dex); weight=rho, region)
            cut = isnothing(at) ? midpoint_layer(mesh,s,region,slice_axis) : at
            plane = slice_plane(mesh, field; axis=slice_axis, at=cut, region)
            export_product(joinpath(stage, "mass_by_density.h5"), hist)
            export_product(joinpath(stage, "slice.h5"), plane)
            resolved = resolve_field(s, field)
            meta = Dict{String,Any}(
                "aegis_version" => string(pkgversion(Aegis)), "source" => s.path, "stream" => chosen,
                "time" => s.time, "time_units" => "code time", "model" => "Newtonian MHD, permeability one",
                "frame" => "simulation frame", "density" => "rho", "velocity" => ["vel1","vel2","vel3"],
                "magnetic_field" => ["Bcc1","Bcc2","Bcc3"], "selection" => "half-open by stored cell center",
                "region" => collect.(collect(region.bounds)), "root_cells" => collect(s.rootsize),
                "selected_cells" => ke.count, "selected_volume" => ke.volume,
                "volume_units" => "code length^3", "energy_units" => "code mass * (code length / code time)^2",
                "kinetic_energy" => ke.integral, "magnetic_energy" => me.integral,
                "workers" => ke.provenance.workers, "chunk_blocks" => chunk_blocks,
                "slice_field" => String(resolved.name), "slice_axis" => slice_axis, "slice_at" => cut,
                "probe_block" => probe.block, "probe_cell" => [first(r) for r in cell],
                "probe_value" => Float64(example_value), "included_mass" => hist.included_weight,
                "excluded_mass" => hist.excluded_weight,
                "histories" => companions.histories, "inputs" => companions.inputs,
                "fields" => [Dict("dataset"=>f.dataset, "slot"=>f.index, "name"=>f.name,
                                  "precision"=>string(f.precision)) for f in s.fields])
            isnothing(s.cycle) || (meta["cycle"] = s.cycle)
            open(io -> TOML.print(io, meta; sorted=true), joinpath(stage,"summary.toml"), "w")
            meta
        end
        mv(stage, dest)
        result
    end
end

function main(args)
    2 <= length(args) <= 3 || error("Usage: julia --project=. examples/mhd_workflow.jl --demo|RUN_DIRECTORY NEW_OUTPUT_DIRECTORY [SLICE_FIELD]")
    output = args[2]
    if args[1] == "--demo"
        result = mktempdir() do dir
            write_demo(dir)
            run_mhd_workflow(dir, output; slice_field=get(args, 3, "tracer"), workers=Threads.nthreads())
        end
        @assert result["kinetic_energy"] == 12600 && result["magnetic_energy"] == 378
        @assert result["included_mass"] ≈ 1008
    else
        result = run_mhd_workflow(args[1], output; slice_field=get(args, 3, nothing), workers=Threads.nthreads())
    end
    println("t=",result["time"],"  cells=",result["selected_cells"],
            "  kinetic energy=",result["kinetic_energy"],"  magnetic energy=",result["magnetic_energy"])
    println("Products: ", abspath(output))
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    MHDWorkflowExample.main(ARGS)
end
