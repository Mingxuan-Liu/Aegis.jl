"""A scalar expression applied independently to native cell samples."""
struct CellMap{F,R,M}
    path::String
    name::String
    units::Union{String,Nothing}
    centering::Symbol
    f::F
    inputs::R
    sources::M
end

"""
    cell_map(f, snapshot, fields...; name, units=nothing)

Define a custom cell-centered scalar, with the given label and optional units.
Each call to `read_field` reads the selected slab of its input fields and
maps the pure scalar function `f` over matching cells, with Float64 arguments
and outputs. Record the formula in `name` or in the calling analysis.
For example, `cell_map(/, s, "H", "C"; name="H/C", units="1")` supplies a
selector for a local ratio. The caller establishes compatibility of input
units and handles the domain of the function.
"""
function cell_map(f, s::Snapshot, fields...; name::AbstractString, units=nothing)
    isempty(fields) && throw(ArgumentError("cell_map requires at least one field"))
    isempty(name) && throw(ArgumentError("give the cell expression a name"))
    inputs = map(x -> _product_ref(s, x), fields)
    sources = map((x,r) -> _product_meta(x,r), fields, inputs)
    CellMap(s.path, String(name), units, :cell, f, inputs, sources)
end

function resolve_field(s::Snapshot, m::CellMap)
    s.path == m.path || throw(ArgumentError("cell expression belongs to another snapshot"))
    foreach(r -> resolve_field(s, r), m.inputs)
    m
end

# Custom maps decide their own domain (for example, a map can handle NaN).
@inline function _cell_map_value(f::F, xs::Vararg{Any,N}) where {F,N}
    Float64(f(map(Float64, xs)...))
end

function read_field(s::Snapshot, m::CellMap, block::Integer;
                    ranges=ntuple(d -> 1:s.blocksize[d], 3))
    resolve_field(s, m)
    arrays = _product_slabs(s, m.inputs, block, ranges)
    _cell_slab(_cell_map_value, m.f, arrays)
end

function _product_ref(s, selector)
    f = resolve_field(s, selector)
    f.centering === :cell || throw(ArgumentError("this product requires cell-centered samples"))
    f
end

# A downstream model can implement the same geometry/reader contract and
# provide ordinary name/unit labels, or specialize this metadata method.
function _product_field_meta(f)
    hasproperty(f,:name) || throw(ArgumentError("a custom product field needs a name for output metadata"))
    units = hasproperty(f,:units) ? something(f.units,"unspecified") : "unspecified"
    (; kind=:user_defined, name=string(f.name), units=string(units))
end

_product_field_meta(f::FieldRef) = (; kind=:stored, f.name, f.dataset, f.index, units="unspecified")
_product_field_meta(f::CellMap) = (; kind=:expression, f.name, units=something(f.units,"unspecified"), inputs=f.sources)

function _product_field_meta(f::NewtonianQuantity)
    eos = f.model.eos
    closure = if isnothing(eos)
        (; kind=:unspecified)
    elseif eos isa IdealGas
        (; kind=:ideal_gas, eos.gamma)
    else
        (; kind=:isothermal, eos.sound_speed)
    end
    (; kind=:newtonian, f.name, f.units, inputs=map(_product_field_meta, f.inputs),
       eos=closure, f.model.velocity_unit_cgs, f.model.mean_particle_mass_g)
end

function _product_meta(selector, ref)
    meta = _product_field_meta(ref)
    selector isa FieldInterpretation || return meta
    merge(meta, (; selector.alias, units=something(selector.units,"unspecified"), selector.description))
end

# At most four selectors per binning pass. Reuse arrays for equal source
# descriptors (e.g. a mass-weighted density PDF) without hashing mesh objects.
function _product_slabs(s, fields::Tuple, block, ranges)
    arrays = Any[]
    for (j,f) in pairs(fields)
        prior = findfirst(i -> fields[i] === f, 1:j-1)
        a = if isnothing(f)
            nothing
        elseif isnothing(prior)
            read_field(s, f, block; ranges)
        else
            arrays[prior]
        end
        push!(arrays, a)
    end
    Tuple(arrays)
end

"""
    BinAxis(field, edges; measure=:linear)

Describe physical, finite, strictly increasing bin edges for one scalar
selector. `:linear` gives density per field unit. `:dex` requires positive
edges and gives density per decade of the same field. Supply edges in
original units in both cases. Finite nonpositive observations on a dex axis
contribute to the excluded weight.
"""
struct BinAxis{F}
    field::F
    edges::Vector{Float64}
    measure::Symbol
    function BinAxis(field, edges; measure::Symbol=:linear)
        e = collect(Float64, edges)
        e isa Vector && length(e) >= 2 && all(isfinite,e) && all(i -> e[i+1]>e[i], 1:length(e)-1) ||
            throw(ArgumentError("bin edges must be a finite, strictly increasing vector"))
        measure in (:linear,:dex) || throw(ArgumentError("bin measure must be :linear or :dex"))
        measure === :dex && e[1] <= 0 && throw(ArgumentError("per-dex bin edges must be positive"))
        new{typeof(field)}(field, e, measure)
    end
end

function _binwidths(a)
    if a.measure === :linear
        widths = diff(a.edges)
    else
        widths = map(1:length(a.edges)-1) do i
            lo, hi = a.edges[i], a.edges[i+1]
            ratio = (hi-lo)/lo
            isfinite(ratio) ? log1p(ratio)/log(10.0) : log10(hi)-log10(lo)
        end
    end
    all(w -> isfinite(w) && w>0, widths) || throw(ArgumentError("unrepresentable bin width"))
    widths
end

function _binindex(edges, x)
    isfinite(x) || throw(DomainError(x, "nonfinite bin coordinate"))
    if x < first(edges) || x > last(edges)
        return 0
    end
    min(searchsortedlast(edges,x), length(edges)-1)
end

function _product_weight(a, i)
    isnothing(a) && return 1.0
    w = Float64(a[i])
    isfinite(w) && w>=0 || throw(DomainError(w,"weight field must be finite and nonnegative"))
    w
end

function _valid_product(a,b)
    t = _weighted_product(a,b)
    isfinite(t) || throw(OverflowError("product exceeds Float64 range"))
    t
end

function _deposit!(totals, corrections, index, value)
    totals[index], corrections[index] = _compensated_add(totals[index], corrections[index], value)
end

# Typed function barrier between HDF5 reads and per-cell arithmetic.
function _bin_block!(totals, tc, nums, nc, axes, coords, samples, weights, block, count_weight)
    rejected = 0.0
    rc = 0.0
    for I in CartesianIndices(first(coords))
        inds = map((a, v) -> _binindex(a.edges,v[I]), axes, coords)
        dv = count_weight ? 1.0 : cell_volume(block,Tuple(I)...)
        isfinite(dv) && dv > 0 || throw(DomainError(dv,"unrepresentable cell volume"))
        w = _valid_product(dv, _product_weight(weights, I))
        if any(iszero, inds)
            rejected, rc = _compensated_add(rejected, rc, w)
            continue
        end
        index = CartesianIndex(inds)
        _deposit!(totals, tc, index, w)
        if !isnothing(nums)
            q = Float64(samples[I])
            isfinite(q) || throw(DomainError(q, "nonfinite conditional value in populated bin"))
            _deposit!(nums, nc, index, _valid_product(w,q))
        end
    end
    rejected + rc
end

function _product_context(s, region)
    (; source_path=s.path, time=s.time, cycle=s.cycle, time_units="code time", coordinate_system=s.coordinates,
       coordinate_units="code length", region_bounds=region.bounds, selection=:half_open_cell_center)
end

# Split mantissas and binary exponents to avoid underflow of a tiny bin
# probability before division by a very narrow bin, or overflow of a joint
# width product. At most two axis widths occur in the supported products.
function _pdf_density(bin_weight, included, widths, I)
    iszero(included) && return NaN
    iszero(bin_weight) && return 0.0
    m,e = frexp(bin_weight)
    m0,e0 = frexp(included)
    m /= m0
    e -= e0
    for d in eachindex(widths)
        mw,ew = frexp(widths[d][I[d]])
        m /= mw
        e -= ew
    end
    result = ldexp(m,e)
    isfinite(result) && !iszero(result) || throw(DomainError(result,"nonzero PDF density exceeds Float64 range"))
    result
end

function _bin_product(mesh, axes, value, weight, region)
    s = mesh.snapshot
    bins = map(a -> _product_ref(s,a.field), axes)
    response = isnothing(value) ? nothing : _product_ref(s,value)
    count_weight = weight === :count
    wf = isnothing(weight) || count_weight ? nothing : _product_ref(s,weight)
    fields = (bins..., response, wf)
    shape = map(a -> length(a.edges)-1, axes)
    widths = map(_binwidths,axes)
    totals, tc = zeros(shape), zeros(shape)
    nums, nc = isnothing(response) ? (nothing,nothing) : (zeros(shape),zeros(shape))
    excluded = 0.0
    exc_c = 0.0
    for b in selectedblocks(mesh; region)
        arrays = _product_slabs(s, fields, b.block, b.ranges)
        coords = ntuple(i -> arrays[i], length(axes))
        rejected = _bin_block!(totals,tc,nums,nc,axes,coords,arrays[end-1],arrays[end],b,count_weight)
        excluded,exc_c = _compensated_add(excluded,exc_c,rejected)
    end
    totals .+= tc
    all(isfinite,totals) || throw(OverflowError("bin total exceeds Float64 range"))
    included = sum(totals)
    excluded += exc_c
    isfinite(included) && isfinite(excluded) || throw(OverflowError("histogram weight exceeds Float64 range"))
    weighting = count_weight ? :count : isnothing(wf) ? :volume : :field_times_volume
    bin_fields = map((a,r) -> _product_meta(a.field,r), axes, bins)
    weight_field = isnothing(wf) ? nothing : _product_meta(weight,wf)
    weight_units = count_weight ? "count" : isnothing(wf) ? "code length^3" : "($(weight_field.units)) * code length^3"
    metadata = (; _product_context(s,region)..., bin_fields, measures=map(a->a.measure,axes),
                  weighting, weight_field, weight_units, last_edge_inclusive=true)
    shared = (; edges=map(a->copy(a.edges),axes), widths, included_weight=included, excluded_weight=excluded)
    if isnothing(response)
        pdf = similar(totals)
        for I in CartesianIndices(totals)
            pdf[I] = _pdf_density(totals[I],included,widths,I)
        end
        axis_units = map((a,f) -> a.measure === :dex ? "dex" : f.units, axes,bin_fields)
        pdf_units = join(("1/($u)" for u in axis_units), " * ")
        return (; kind=:distribution, shared..., totals, pdf,
                   metadata=merge(metadata,(; normalization=:included_weight, pdf_units)))
    end
    nums .+= nc
    all(isfinite,nums) || throw(OverflowError("conditional numerator exceeds Float64 range"))
    means = map((n,d) -> iszero(d) ? NaN : n/d, nums, totals)
    all(i -> iszero(totals[i]) || isfinite(means[i]), eachindex(nums)) || throw(OverflowError("conditional mean exceeds Float64 range"))
    value_field = _product_meta(value,response)
    units = (; numerator="($weight_units) * ($(value_field.units))", denominator=weight_units, mean=value_field.units)
    (; kind=:conditional_stats, shared..., numerator=nums, denominator=totals, mean=means,
       metadata=merge(metadata, (; value_field, normalization=:per_bin_weight, units)))
end

"""
    distribution(mesh, xaxis[, yaxis]; weight=nothing, region=Region())

Return Float64 bin `totals` and a `pdf` normalized by the included weight and
each axis's bin width. A cell has weight dV by default, 1 for `:count`, or
weight_field × dV for a supplied selector. Each selected cell contributes
once. All bins are lower-inclusive; the last also includes its upper edge.
Empty bins have zero density when some weight is included; an entirely
zero-weight domain has undefined (NaN) density. `excluded_weight` counts
weight outside any axis once. Axes follow argument order in Julia. Units of
the PDF are inverse products of axis units or decades. See docs/products.md.
"""
distribution(mesh::CartesianMesh, x::BinAxis; weight=nothing, region::Region=Region()) =
    _bin_product(mesh,(x,),nothing,weight,region)
distribution(mesh::CartesianMesh, x::BinAxis, y::BinAxis; weight=nothing, region::Region=Region()) =
    _bin_product(mesh,(x,y),nothing,weight,region)

"""
    conditional_stats(mesh, axis, value; weight=nothing, region=Region())

Bin a scalar value by another scalar using the edge and weighting rules of
`distribution`. Return `numerator = Σ w value`, `denominator = Σ w`, and
`mean` per bin. Empty or zero-weight bins have zero sums and NaN means.
A nonfinite response inside the bin domain raises an error. The numerator
has units of weight times response, and the mean has response units.
"""
conditional_stats(mesh::CartesianMesh, a::BinAxis, value; weight=nothing, region::Region=Region()) =
    _bin_product(mesh,(a,),value,weight,region)

function _global_axis(mesh, d)
    segments = mesh.axes[d]
    centers = Vector{Float64}(undef,mesh.snapshot.rootsize[d])
    faces = Vector{Float64}(undef,length(centers)+1)
    offset = 1
    for segment in segments
        n = length(segment.centers)
        centers[offset:offset+n-1] .= segment.centers
        faces[offset:offset+n] .= segment.faces
        offset += n
    end
    (; centers, faces)
end

function _plane_grid(mesh, region, axis)
    axis in 1:3 || throw(ArgumentError("axis must be 1, 2 or 3"))
    all_axes = ntuple(d -> _global_axis(mesh,d),3)
    keep = ntuple(3) do d
        lo,hi = region.bounds[d]
        centers = all_axes[d].centers
        searchsortedfirst(centers,lo):searchsortedfirst(centers,hi)-1
    end
    any(isempty, keep) && throw(ArgumentError("region has no native cells"))
    transverse = Tuple(d for d in 1:3 if d != axis)
    centers = map(d -> all_axes[d].centers[keep[d]],transverse)
    faces = map(d -> all_axes[d].faces[first(keep[d]):last(keep[d])+1],transverse)
    (; all_axes, keep, transverse, centers, faces)
end

function _pixel_ranges(mesh,b,grid)
    loc = mesh.snapshot.locations
    bs = mesh.snapshot.blocksize
    map(grid.transverse) do d
        offset = loc[d,b.block]*bs[d] - first(grid.keep[d]) + 1
        (offset+first(b.ranges[d])):(offset+last(b.ranges[d]))
    end
end

"""
    slice_plane(mesh, field; axis, at, region=Region())

Read one native cell layer intersected by a coordinate plane. On an internal
face, select the cell on the positive side. The outermost upper face lies
outside the half-open domain. Selection uses stored centers and the chosen
layer must belong to the region. Output `values` has field units, and its
two axes are the surviving physical axes in ascending order. `faces` and
`centers` locate each pixel in code lengths. No spatial interpolation is used.
"""
function slice_plane(mesh::CartesianMesh, field; axis::Int, at::Real, region::Region=Region())
    s = mesh.snapshot
    ref = _product_ref(s,field)
    grid = _plane_grid(mesh,region,axis)
    plane = Float64(at)
    isfinite(plane) || throw(ArgumentError("slice coordinate must be finite"))
    layer = searchsortedlast(grid.all_axes[axis].faces,plane)
    layer in grid.keep[axis] || throw(ArgumentError("plane does not intersect a selected cell layer"))
    values = Array{Float64}(undef,map(length,grid.centers))
    for b in selectedblocks(mesh;region)
        i = layer - s.locations[axis,b.block]*s.blocksize[axis]
        i in b.ranges[axis] || continue
        slab = read_field(s,ref,b.block; ranges=Base.setindex(b.ranges,i:i,axis))
        all(isfinite,slab) || throw(DomainError(field,"nonfinite field in selected slice"))
        pixels = _pixel_ranges(mesh,b,grid)
        values[pixels...] .= dropdims(slab;dims=axis)
    end
    layer_faces = grid.all_axes[axis].faces[layer:layer+1]
    metadata = (; _product_context(s,region)..., field=_product_meta(field,ref), axis,
                  transverse_axes=grid.transverse, at=plane, layer, layer_faces,
                  layer_center=grid.all_axes[axis].centers[layer], method=:piecewise_constant_slice)
    (; kind=:slice, values, grid.faces, grid.centers, metadata)
end

function _project_block!(nums,nc,den,dc,values,weights,widths,axis,transverse,pixels)
    for I in CartesianIndices(values)
        q = Float64(values[I])
        isfinite(q) || throw(DomainError(q,"nonfinite field in projection"))
        dl = widths[I[axis]]
        w = _valid_product(dl, _product_weight(weights,I))
        index = CartesianIndex(pixels[1][I[transverse[1]]],pixels[2][I[transverse[2]]])
        _deposit!(nums,nc,index,_valid_product(q,w))
        isnothing(den) || _deposit!(den,dc,index,w)
    end
end

"""
    project(mesh, field; axis, mode=:integral, weight=nothing, region=Region())

Compute a native Cartesian line integral (`Σ q Δl`) or a mean along each
ray (`mode=:mean`, `Σ w q Δl / Σ w Δl`). The unweighted mean uses w=1.
A supplied nonnegative weight is allowed for the mean. Zero weight on a ray
gives NaN. A pixel covers one native transverse cell with the two recorded
face intervals. Multiply the line integral by pixel area before summing
to recover a volume integral. Result `values` has field × length units for
integrals or field units for means. `path_length` is the common ray length
in code units. `denominator` is a weight × length map for weighted means,
the scalar path length for unweighted means, or nothing for integrals.
"""
function project(mesh::CartesianMesh, field; axis::Int, mode::Symbol=:integral, weight=nothing,
                 region::Region=Region())
    mode in (:integral,:mean) || throw(ArgumentError("projection mode must be :integral or :mean"))
    mode === :integral && !isnothing(weight) && throw(ArgumentError("a weight field applies to mode=:mean"))
    s = mesh.snapshot
    ref = _product_ref(s,field)
    wf = isnothing(weight) ? nothing : _product_ref(s,weight)
    grid = _plane_grid(mesh,region,axis)
    shape = map(length, grid.centers)
    nums,nc = zeros(shape),zeros(shape)
    den,dc = isnothing(wf) ? (nothing,nothing) : (zeros(shape),zeros(shape))
    for b in selectedblocks(mesh;region)
        arrays = _product_slabs(s,(ref,wf),b.block,b.ranges)
        pixels = _pixel_ranges(mesh,b,grid)
        _project_block!(nums,nc,den,dc,arrays...,b.widths[axis],axis,grid.transverse,pixels)
    end
    nums .+= nc
    all(isfinite,nums) || throw(OverflowError("line integral exceeds Float64 range"))
    faces = grid.all_axes[axis].faces
    r = grid.keep[axis]
    path = faces[last(r)+1] - faces[first(r)]
    isfinite(path) && path>0 || throw(OverflowError("unrepresentable path length"))
    if !isnothing(den)
        den .+= dc
        all(isfinite,den) || throw(OverflowError("weight integral exceeds Float64 range"))
        values = map((n,d) -> iszero(d) ? NaN : n/d, nums,den)
    elseif mode === :mean
        den = path
        values = nums ./ den
    else
        values = nums
    end
    any(isinf,values) && throw(OverflowError("projection exceeds Float64 range"))
    source = _product_meta(field,ref)
    weight_field = isnothing(wf) ? nothing : _product_meta(weight,wf)
    values_units = mode === :integral ? "($(source.units)) * code length" : source.units
    denominator_units = isnothing(wf) ? "code length" : "($(weight_field.units)) * code length"
    metadata = (; _product_context(s,region)..., field=source, axis, values_units, denominator_units,
                  transverse_axes=grid.transverse, method=mode, ray_faces=(faces[first(r)],faces[last(r)+1]),
                  weighting=isnothing(wf) ? :path_length : :field_times_path_length, weight_field)
    (; kind=:projection, values, grid.faces, grid.centers, path_length=path, denominator=den, metadata)
end

# HDF5 named subgroups preserve provenance without display-dependent repr
# parsing. An absent optional property is omitted throughout schema version 1.
function _write_product_node(parent,key,value)
    if isnothing(value)
        return
    elseif value isa NamedTuple || value isa Tuple
        group = HDF5.create_group(parent,key)
        try
            for (k,v) in pairs(value)
                _write_product_node(group,string(k),v)
            end
        finally
            close(group)
        end
    else
        parent[key] = value isa Symbol ? String(value) : value
    end
end

"""
    export_product(path, product)

Write a distribution, conditional statistic, slice or projection to a new
HDF5 file. Root datasets/groups mirror the result's named properties. Tuples
use one-based numbered child names. Optional absent metadata are omitted.
The root attributes record the format and Aegis version. Array axes in
Julia follow the product's order; a C/Python HDF5 reader sees these array
dimensions reversed. Use the exported axis edges, coordinates, units,
weighting and selection to interpret the numerical result. Existing outputs
are protected. See docs/products.md for the schema.
"""
function export_product(path::AbstractString, product::NamedTuple)
    hasproperty(product,:kind) && product.kind in (:distribution,:conditional_stats,:slice,:projection) ||
        throw(ArgumentError("expected an Aegis numerical product"))
    dest = abspath(path)
    ispath(dest) && throw(ArgumentError("product path already exists: $dest"))
    mktemp(dirname(dest)) do tmp, io
        close(io)
        HDF5.h5open(tmp,"w") do h
            HDF5.attributes(h)["aegis_product_schema"] = 1
            HDF5.attributes(h)["aegis_version"] = string(Base.pkgversion(@__MODULE__))
            for (key,value) in pairs(product)
                _write_product_node(h,string(key),value)
            end
        end
        # Atomic creation on a single filesystem, including simultaneous
        # attempts by two writers. mktemp removes its temporary link on exit.
        hardlink(tmp,dest)
    end
    return dest
end
