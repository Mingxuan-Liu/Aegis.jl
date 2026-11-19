"""
    FieldRef

A field's identity in one snapshot: path, original dataset and name, and
one-based variable index within the dataset. `precision` is the element type
on disk. `centering` is `:cell` for the supported Athena++ cell-output
contract, or `:unknown` when the writer can also pack face samples under
the same label. Units and physical definitions are separate interpretations.
"""
struct FieldRef
    path::String
    dataset::String
    index::Int
    name::String
    precision::DataType
    centering::Symbol
end

"""
    Snapshot

Native mesh metadata and read-only dataset handles owned by `open_snapshot`.
Sizes follow physical axis order (x1,x2,x3); `locations[:,b]` are the stored
zero-based logical coordinates at `levels[b]`. `rootgrid[d]` retains the
minimum, maximum and spacing-ratio attribute for axis `d`. `attributes` holds
the original root attributes in their decoded types. The reported block size
includes any ghost or collapsed zones present in the file.
"""
struct Snapshot
    path::String
    file::HDF5.File
    attributes::Dict{String,Any}
    time::Float64
    cycle::Union{Nothing,Int}
    coordinates::String
    nblocks::Int
    blocksize::NTuple{3,Int}
    rootsize::NTuple{3,Int}
    rootgrid::NTuple{3,NTuple{3,Float64}}
    maxlevel::Int
    levels::Vector{Int}
    locations::Matrix{Int}
    fields::Vector{FieldRef}
    datasets::Dict{String,HDF5.Dataset}
    faces::NTuple{3,HDF5.Dataset}
    centers::NTuple{3,HDF5.Dataset}
    identities::Dict{Tuple{String,Int},FieldRef}
    names::Dict{String,Vector{FieldRef}}
end

function _required(attrs, key)
    haskey(attrs, key) || throw(ArgumentError("missing required .athdf attribute '$key'"))
    attrs[key]
end

function _triple(attrs, key, ::Type{T}) where T
    v = _required(attrs, key)
    v isa AbstractVector && length(v) == 3 && eltype(v) <: T ||
        throw(ArgumentError("'$key' must have three $T entries in physical axis order"))
    Tuple(v)
end

function _integer(value, key)
    value isa Integer || throw(ArgumentError("'$key' must be a scalar integer"))
    return Int(value)
end

function _dataset(file, key, shape)
    haskey(file, key) || throw(ArgumentError("missing .athdf dataset '$key'"))
    ds = file[key]
    ds isa HDF5.Dataset || throw(ArgumentError("'$key' is not a dataset"))
    size(ds) == shape || throw(ArgumentError(
        "dataset '$key' has Julia dimensions $(size(ds)); expected $shape from .athdf metadata"))
    return ds
end

function _coordinates_dataset(file, name, shape)
    ds = _dataset(file, name, shape)
    eltype(ds) in (Float32, Float64) || throw(ArgumentError(
        "coordinate dataset '$name' requires floating-point code lengths"))
    return ds
end

function _make_snapshot(file, path)
    attrs = Dict{String,Any}(key => HDF5.read_attribute(file, key)
                             for key in keys(HDF5.attributes(file)))
    t = _required(attrs, "Time")
    t isa Real && isfinite(t) || throw(ArgumentError("'Time' must be a finite scalar in code units"))
    cycle = haskey(attrs, "NumCycles") ? _integer(attrs["NumCycles"], "NumCycles") : nothing
    coord = _required(attrs, "Coordinates")
    coord isa AbstractString && !isempty(coord) ||
        throw(ArgumentError("'Coordinates' must specify a coordinate system"))
    nblocks = _integer(_required(attrs, "NumMeshBlocks"), "NumMeshBlocks")
    nblocks > 0 || throw(ArgumentError("'NumMeshBlocks' must be positive"))
    blocksize = Int.(_triple(attrs, "MeshBlockSize", Integer))
    rootsize = Int.(_triple(attrs, "RootGridSize", Integer))
    all(>(0), blocksize) && all(>(0), rootsize) ||
        throw(ArgumentError("mesh block and root-grid sizes must be positive"))
    rootgrid = ntuple(d -> Float64.(_triple(attrs, "RootGridX$d", Real)), 3)
    maxlevel = _integer(_required(attrs, "MaxLevel"), "MaxLevel")
    levels_ds = _dataset(file, "Levels", (nblocks,))
    locations_ds = _dataset(file, "LogicalLocations", (3, nblocks))
    eltype(levels_ds) <: Integer && eltype(locations_ds) <: Integer ||
        throw(ArgumentError("Levels and LogicalLocations must contain integers"))
    levels = Int.(read(levels_ds))
    locations = Int.(read(locations_ds))
    close(levels_ds)
    close(locations_ds)
    all(l -> 0 <= l <= maxlevel, levels) && all(>=(0), locations) ||
        throw(ArgumentError("invalid levels or logical locations in .athdf mesh metadata"))
    faces = ntuple(d -> _coordinates_dataset(file, "x$(d)f", (blocksize[d]+1, nblocks)), 3)
    centers = ntuple(d -> _coordinates_dataset(file, "x$(d)v", (blocksize[d], nblocks)), 3)

    dataset_names = _required(attrs, "DatasetNames")
    counts = _required(attrs, "NumVariables")
    varnames = _required(attrs, "VariableNames")
    dataset_names isa AbstractVector{<:AbstractString} && !isempty(dataset_names) &&
        allunique(dataset_names) && all(!isempty, dataset_names) ||
        throw(ArgumentError("'DatasetNames' must list distinct, nonempty dataset names"))
    counts isa AbstractVector{<:Integer} && length(counts) == length(dataset_names) &&
        all(>(0), counts) || throw(ArgumentError("'NumVariables' requires a positive count per dataset"))
    varnames isa AbstractVector{<:AbstractString} && sum(Int, counts) == length(varnames) ||
        throw(ArgumentError("'VariableNames' length disagrees with 'NumVariables'"))

    datasets = Dict{String,HDF5.Dataset}()
    fields = FieldRef[]
    identities = Dict{Tuple{String,Int},FieldRef}()
    names = Dict{String,Vector{FieldRef}}()
    offset = 0
    for (dataset_name, count) in zip(dataset_names, counts)
        ds = _dataset(file, dataset_name, (blocksize..., nblocks, Int(count)))
        precision = eltype(ds)
        precision in (Float32, Float64) || throw(ArgumentError(
            "unsupported element type $precision for '$dataset_name'; expected Float32 or Float64"))
        datasets[dataset_name] = ds
        for i in 1:Int(count)
            name = varnames[offset+i]
            # outputs.cpp can pack the left face of each cell into the same
            # shape as cell samples, even when XDMF declares Center="Cell".
            centering = name in ("B1", "B2", "B3") ? :unknown : :cell
            ref = FieldRef(path, dataset_name, i, name, precision, centering)
            push!(fields, ref)
            identities[(dataset_name, i)] = ref
            push!(get!(Vector{FieldRef}, names, name), ref)
        end
        offset += count
    end
    return Snapshot(path, file, attrs, Float64(t), cycle, coord, nblocks, blocksize, rootsize,
                    rootgrid, maxlevel, levels, locations, fields, datasets, faces, centers,
                    identities, names)
end

"""
    open_snapshot(f, path)

Open an Athena++ snapshot read-only for the duration of `f(snapshot)`. Read
root attributes, field descriptors, levels and logical locations; field data
and coordinate arrays remain on disk until requested. All HDF5 handles are
closed when the callback returns or throws. Arrays returned by reads are owned
by the caller. Use this handle from one I/O task at a time.

Requires conventional Athena++ HDF5 field groups with C storage order
(variable, block, x3, x2, x1). Fields may be omitted, renamed, repeated, or
spread over datasets. Unscaled Float32 and Float64 fields are supported.
See the Athena++ writer `src/outputs/athena_hdf5.cpp`, `WriteOutputFile`,
and the reference reader `vis/python/athena_read.py`, `athdf(raw=True)`.
"""
function open_snapshot(f::Function, path::AbstractString)
    fullpath = abspath(path)
    HDF5.h5open(fullpath, "r") do file
        f(_make_snapshot(file, fullpath))
    end
end

"""
    resolve_field(snapshot, selector)

Resolve a unique original name, a source key `(dataset, one_based_index)`,
or a descriptor obtained from the same snapshot path. Ambiguous names require
the source key so two custom fields can always be distinguished.
"""
function resolve_field(s::Snapshot, name::AbstractString)
    matches = get(s.names, name, nothing)
    isnothing(matches) && throw(ArgumentError("field '$name' absent in $(s.path); inspect snapshot.fields"))
    length(matches) == 1 && return only(matches)
    candidates = [(ref.dataset, ref.index) for ref in matches]
    throw(ArgumentError("field '$name' is ambiguous; select one of $candidates"))
end

function resolve_field(s::Snapshot, key::Tuple{AbstractString,Integer})
    return get(s.identities, key) do
        throw(ArgumentError("field identity $key absent in $(s.path); inspect snapshot.fields"))
    end
end

function resolve_field(s::Snapshot, ref::FieldRef)
    ref.path == s.path || throw(ArgumentError("field descriptor belongs to a different snapshot path"))
    found = resolve_field(s, (ref.dataset, ref.index))
    found == ref || throw(ArgumentError("field descriptor no longer matches snapshot metadata"))
    return found
end

function _checkblock(s, b)
    isopen(s.file) || throw(ArgumentError("snapshot handle closed; read inside open_snapshot's do block"))
    1 <= b <= s.nblocks || throw(ArgumentError("block $b outside 1:$(s.nblocks)"))
end

function _slab_indices(s::Snapshot, ranges)
    ranges isa Tuple && length(ranges) == 3 ||
        throw(ArgumentError("ranges must be a tuple of three contiguous cell-index ranges"))
    return ntuple(3) do d
        r = ranges[d]
        r isa AbstractUnitRange{<:Integer} && !isempty(r) && first(r) >= 1 && last(r) <= s.blocksize[d] ||
            throw(ArgumentError("range $r must select cells within 1:$(s.blocksize[d]) on axis $d"))
        # HDF5.jl hyperslabs use native Int indices; bounds were checked before narrowing.
        Int(first(r)):Int(last(r))
    end
end

"""
    read_field(snapshot, selector, block; ranges=...)

Read one native block (one-based storage index) or one contiguous hyperslab
in each physical axis. `ranges` is a tuple of three nonempty, one-based unit
ranges in (x1,x2,x3) order. The returned array has those same axes, including
singleton axes, and retains the dataset precision. No units, ghost trimming,
interpolation, or refinement restriction are applied.
"""
function read_field(s::Snapshot, selector, block::Integer;
                    ranges=ntuple(d -> 1:s.blocksize[d], 3))
    _checkblock(s, block)
    ref = resolve_field(s, selector)
    indices = _slab_indices(s, ranges)
    # HDF5.jl reverses C axis order to accommodate Julia's column-major layout.
    return s.datasets[ref.dataset][indices..., Int(block), ref.index]
end

"""
    block_coordinates(snapshot, block)

Read the stored `faces` and `centers` coordinate vectors for one native block
without reconstructing a global mesh. Both tuples are in physical (x1,x2,x3)
order and retain source precision. For Cartesian coordinates all three axes are
in code length units; other coordinate systems retain their native convention.
"""
function block_coordinates(s::Snapshot, block::Integer)
    _checkblock(s, block)
    return (faces=ntuple(d -> s.faces[d][:, Int(block)], 3),
            centers=ntuple(d -> s.centers[d][:, Int(block)], 3))
end

"""One snapshot's absolute path, validated code time and optional solver cycle."""
struct SnapshotRecord
    path::String
    time::Float64
    cycle::Union{Nothing,Int}
end

"""
    discover_snapshots(directory)

Inspect each immediate `.athdf` file and return a dictionary of stream keys
to snapshot records. For `prefix.00042.athdf`, the key is `prefix`; other
basenames are individual keys including the extension. Each stream is sorted
by stored time, then path. Equal times remain separate records. No cross-stream
synchronization or restart resolution is implied. A filename that collides
with another stream's key raises an ambiguity error. Files with invalid snapshot
metadata raise an error identifying the path. No mesh field data is loaded.
"""
function discover_snapshots(directory::AbstractString)
    streams = Dict{String,Vector{SnapshotRecord}}()
    numbered_stream = Dict{String,Bool}()
    for path in readdir(directory; join=true)
        isfile(path) && endswith(path, ".athdf") || continue
        name = basename(path)
        numbered = match(r"^(.+)\.[0-9]+\.athdf$", name)
        key = isnothing(numbered) ? name : String(numbered.captures[1])
        get!(numbered_stream, key, !isnothing(numbered)) == !isnothing(numbered) ||
            throw(ArgumentError("unnumbered snapshot and numbered stream share key '$key'"))
        record = try
            open_snapshot(path) do s
                SnapshotRecord(s.path, s.time, s.cycle)
            end
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("cannot catalog '$path': $(sprint(showerror, err))"))
        end
        push!(get!(Vector{SnapshotRecord}, streams, key), record)
    end
    for records in values(streams)
        sort!(records; by=r -> (r.time, r.path))
    end
    return streams
end
