using HDF5

# Manufactured native blocks. Deliberately independent of the Aegis reader.
# Writer order in HDF5.jl is (x1, x2, x3, block, variable), reversed on disk.
function write_fixture(path; order=[2, 1], fields=("marker", "rho", "tracer"),
                       blocksize=(3, 2, 2), time=1.25, cycle=8)
    nx, ny, nz = blocksize
    nb = length(order)
    sort(order) == 1:nb || error("order must permute 1:nb")
    rootsize = (nb * nx, ny, nz)
    groups = filter(pair -> !isempty(last(pair)), [
        "gas" => [s for s in fields if s == "rho"],
        "extras" => [s for s in fields if s != "rho"],
    ])
    h5open(path, "w") do h
        a = attributes(h)
        a["Coordinates"] = "cartesian"
        a["Time"] = Float64(time)
        a["NumCycles"] = Int64(cycle)
        a["NumMeshBlocks"] = Int32(nb)
        a["MaxLevel"] = Int32(0)
        a["RootGridSize"] = Int32[rootsize...]
        a["MeshBlockSize"] = Int32[blocksize...]
        a["DatasetNames"] = first.(groups)
        a["NumVariables"] = Int32[length(names) for (_, names) in groups]
        a["VariableNames"] = reduce(vcat, last.(groups))

        h["Levels"] = zeros(Int32, nb)
        logical = zeros(Int64, 3, nb)
        logical[1, :] = order .- 1
        h["LogicalLocations"] = logical
        for d in 1:3
            # Integer cell widths 1,2,... make weighted answers hand-calculable.
            faces = Float64[i * (i + 1) ÷ 2 for i in 0:rootsize[d]]
            a["RootGridX$d"] = [first(faces), last(faces), -1.0] # user mesh generator
            n = blocksize[d]
            blockfaces = Matrix{Float64}(undef, n + 1, nb)
            for (b, loc) in enumerate(order)
                offset = d == 1 ? (loc - 1) * nx : 0
                blockfaces[:, b] = faces[offset+1:offset+n+1]
            end
            h["x$(d)f"] = blockfaces
            h["x$(d)v"] = (blockfaces[1:end-1, :] + blockfaces[2:end, :]) / 2
        end
        for (dataset, names) in groups
            data = Array{Float32}(undef, nx, ny, nz, nb, length(names))
            for (f, name) in enumerate(names), (b, loc) in enumerate(order)
                for k in 1:nz, j in 1:ny, i in 1:nx
                    gi = (loc - 1) * nx + i
                    value = if name == "marker"
                        gi + 10j + 100k
                    elseif name == "rho"
                        1 + gi
                    else
                        gi
                    end
                    data[i, j, k, b, f] = value
                end
            end
            h[dataset] = data
        end
    end
    path
end

"""Write a small, synthetic run with a passive scalar on two Cartesian blocks."""
function example_run(directory)
    input = joinpath(directory, "athinput.tracer")
    history = joinpath(directory, "demo.hst")
    snapshot = joinpath(directory, "demo.out2.00001.athdf")
    write(input, """
    <job>
    problem_id = demo
    <problem>
    tracer_scale = 2.5   # user-defined passive scalar parameter
    <par_end>
    """)
    write(history, """
    # Athena++ history data
    # [1]=time [2]=dt [3]=tracer integral
    1.25 0.01 819
    """)
    write_fixture(snapshot)
    (; input, history, snapshot)
end
