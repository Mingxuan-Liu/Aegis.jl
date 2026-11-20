using HDF5

# Pack manufactured full-grid samples into an Athena++-layout file. Block
# storage order is deliberately reversed relative to logical position.
function write_uniform_snapshot(path, fields; lengths, tiles=(1,1,1),
                                names=["component$i" for i in eachindex(fields)])
    dims = size(first(fields))
    bs = dims .÷ tiles
    @assert bs .* tiles == dims
    logical = reverse(vec(collect(CartesianIndices(tiles))))
    nb = length(logical)
    @assert length(names) == length(fields) && all(f -> size(f) == dims, fields)
    h5open(path, "w") do h
        a = attributes(h)
        a["Coordinates"] = "cartesian"
        a["Time"] = 2.0
        a["NumMeshBlocks"] = nb
        a["MaxLevel"] = 0
        a["RootGridSize"] = collect(dims)
        a["MeshBlockSize"] = collect(bs)
        a["DatasetNames"] = ["scalars"]
        a["NumVariables"] = [length(fields)]
        a["VariableNames"] = names
        h["Levels"] = zeros(Int32, nb)
        h["LogicalLocations"] = [c[d]-1 for d in 1:3, c in logical]
        for d in 1:3
            lo = -lengths[d]/2
            dx = lengths[d]/dims[d]
            a["RootGridX$d"] = [lo, -lo, 1.0]
            f = [lo + dx*((c[d]-1)*bs[d] + i) for i in 0:bs[d], c in logical]
            h["x$(d)f"] = f
            h["x$(d)v"] = (f[1:end-1,:] + f[2:end,:])/2
        end
        data = Array{Float64}(undef, bs..., nb, length(fields))
        for (b,c) in enumerate(logical)
            ranges = ntuple(d -> ((c[d]-1)*bs[d]+1):(c[d]*bs[d]), 3)
            for j in eachindex(fields)
                data[:,:,:,b,j] .= view(fields[j], ranges...)
            end
        end
        h["scalars"] = data
    end
    path
end
