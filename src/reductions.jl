"""
    Region(; x1=(-Inf, Inf), x2=(-Inf, Inf), x3=(-Inf, Inf))

Select cells whose stored centers lie in each half-open coordinate interval
`lo ≤ x < hi`, in simulation length units. A selected cell contributes its
whole volume, including any part outside the interval. Equal bounds select
nothing. The bounds do not wrap at periodic boundaries.
"""
struct Region
    bounds::NTuple{3,Tuple{Float64,Float64}}
    function Region(; x1=(-Inf, Inf), x2=(-Inf, Inf), x3=(-Inf, Inf))
        bounds = map((x1, x2, x3)) do interval
            length(interval) == 2 || throw(ArgumentError("each region axis needs two bounds"))
            lo, hi = Float64.(interval)
            lo <= hi || throw(ArgumentError("region bounds must be ordered and cannot be NaN"))
            (lo, hi)
        end
        new(bounds)
    end
end

# All transverse copies of a Cartesian coordinate segment must agree. Keep one
# segment per logical position along an axis instead of one copy per block.
struct AxisSegment
    faces::Vector{Float64}
    centers::Vector{Float64}
    widths::Vector{Float64}
end

"""
    CartesianMesh(snapshot)

Verify and retain the coordinates of a complete level-zero Cartesian grid for
repeated selections within the snapshot's open lifetime. Accept nonuniform
spacing and physical singleton thicknesses. Require the full logical block
lattice, matching face coordinates at every seam and root boundary, and
strictly increasing faces with a center inside every cell. Use this geometry
only when saved Cartesian lengths define the desired measure. A transformed
or expanding-box physical metric requires a separate conversion. Refinement, ghost
layers, slices and summed outputs that do not tile the root grid are rejected.

Cell volume follows Athena++ `Coordinates::GetCellVolume` for Cartesian
geometry: `Δx1 Δx2 Δx3` from output faces. See Stone et al. (2020), §2.1–2.2,
doi:10.3847/1538-4365/ab929b. Geometrical information retained here scales
with the three root axis lengths; coverage verification also visits every
block and uses one bit per logical location.
"""
struct CartesianMesh{S}
    snapshot::S
    axes::NTuple{3,Vector{AxisSegment}}
end

function CartesianMesh(s::Snapshot)
    s.coordinates == "cartesian" ||
        throw(ArgumentError("native-cell reductions require cartesian coordinates, got $(s.coordinates)"))
    s.maxlevel == 0 && all(iszero, s.levels) ||
        throw(ArgumentError("native-cell reductions currently require level-zero blocks; use raw block access for refinement"))
    all(d -> s.rootsize[d] % s.blocksize[d] == 0, 1:3) ||
        throw(ArgumentError("stored cells do not tile RootGridSize (possible ghost, sliced, or summed output)"))
    n = s.rootsize .÷ s.blocksize
    prod(big.(n)) == s.nblocks ||
        throw(ArgumentError("stored blocks do not cover the root grid exactly once; check for ghosts, slices, or missing blocks"))
    seen = falses(n)
    segments = ntuple(d -> Vector{AxisSegment}(undef, n[d]), 3)
    for b in 1:s.nblocks
        loc = ntuple(d -> s.locations[d,b] + 1, 3)
        all(d -> 1 <= loc[d] <= n[d], 1:3) ||
            throw(ArgumentError("block $b has an out-of-grid logical location"))
        seen[loc...] && throw(ArgumentError("duplicate logical location at block $b"))
        seen[loc...] = true
        coords = block_coordinates(s, b)
        for d in 1:3
            fs, cs = coords.faces[d], coords.centers[d]
            t = loc[d]
            if isassigned(segments[d], t)
                ref = segments[d][t]
                fs == ref.faces && cs == ref.centers ||
                    throw(ArgumentError("inconsistent Cartesian coordinates on axis $d of block $b"))
            else
                length(fs) == s.blocksize[d] + 1 && length(cs) == s.blocksize[d] ||
                    throw(ArgumentError("face or center count disagrees with block size on axis $d"))
                # The writer can store root metadata in double and coordinates
                # in single precision. Compare after the writer's final rounding.
                if (t == 1 && first(fs) != convert(eltype(fs), s.rootgrid[d][1])) ||
                   (t == n[d] && last(fs) != convert(eltype(fs), s.rootgrid[d][2]))
                    throw(ArgumentError("outer faces of block $b do not match axis $d root bounds"))
                end
                faces = Float64.(fs)
                widths = diff(faces)
                all(w -> isfinite(w) && w > 0, widths) &&
                    all(i -> faces[i] < cs[i] < faces[i+1], eachindex(cs)) ||
                    throw(ArgumentError("invalid face spacing or cell centers on axis $d of block $b"))
                segments[d][t] = AxisSegment(faces, Float64.(cs), widths)
            end
        end
    end
    for d in 1:3, t in 2:n[d]
        last(segments[d][t-1].faces) == first(segments[d][t].faces) ||
            throw(ArgumentError("gap or overlap between Cartesian blocks on axis $d"))
    end
    CartesianMesh(s, segments)
end

struct SelectedBlock{W}
    block::Int
    ranges::NTuple{3,UnitRange{Int}}
    widths::W
end

function _select_block(mesh, b, region)
    axes = ntuple(d -> mesh.axes[d][mesh.snapshot.locations[d,b] + 1], 3)
    ranges = ntuple(3) do d
        lo, hi = region.bounds[d]
        centers = axes[d].centers
        searchsortedfirst(centers, lo):(searchsortedfirst(centers, hi) - 1)
    end
    widths = ntuple(d -> view(axes[d].widths, ranges[d]), 3)
    SelectedBlock(b, ranges, widths)
end

"""
    selectedblocks(mesh; region=Region())

Iterate nonempty native blocks in source-file order. Each selection supplies
`block` (one-based storage index), `ranges` (three one-based index ranges
within that native block), and `widths` (shared views of the selected cell
widths in x1,x2,x3 order). An empty region yields no blocks. Pass the ranges
to `read_field` for a matching hyperslab. Treat coordinate views as immutable.
"""
function selectedblocks(mesh::CartesianMesh; region::Region=Region())
    blocks = (_select_block(mesh, b, region) for b in 1:mesh.snapshot.nblocks)
    Iterators.filter(b -> all(r -> !isempty(r), b.ranges), blocks)
end

"""
    cell_volume(selection, i, j, k)

Cartesian cell volume in code length cubed from the three stored face
separations, including the stored thickness of singleton axes. Indices are
one-based offsets **within the selection** returned by `selectedblocks`.
Volumes inferred from rounded output coordinates can differ from the solver's
internal precision.
"""
@inline function cell_volume(b::SelectedBlock, i::Int, j::Int, k::Int)
    b.widths[1][i] * b.widths[2][j] * b.widths[3][k]
end

# A shifted mean retains small fluctuations on large constant backgrounds.
# norm2 is the *unnormalized* weighted central second moment. The shifted mean
# and norm2 are merged by the parallel weighted-variance identity.
mutable struct CellTotals
    count::Int
    volume::Float64
    volume_correction::Float64
    integral::Float64
    correction::Float64
    numerator::Float64
    numerator_correction::Float64
    minimum::Float64
    maximum::Float64
    weight::Float64
    weight_correction::Float64
    origin::Float64
    offset::Float64
    norm2::Float64
end

CellTotals() = CellTotals(0, 0, 0, 0, 0, 0, 0, Inf, -Inf, 0, 0, 0, 0, 0)

# Compensated summation matters for signed fields with strong cancellation.
@inline function _compensated_add(value, correction, term)
    s = value + term
    residual = abs(value) >= abs(term) ? (value - s) + term : (term - s) + value
    return (s, correction + residual)
end

@inline function _add_moments!(a, w, origin, offset, norm2, wcor=0.0)
    w == 0 && return
    if a.weight == 0
        a.weight, a.weight_correction = w, wcor
        a.origin, a.offset, a.norm2 = origin, offset, norm2
        return
    end
    oldweight = a.weight + a.weight_correction
    a.weight, a.weight_correction = _compensated_add(a.weight, a.weight_correction, w)
    a.weight, a.weight_correction = _compensated_add(a.weight, a.weight_correction, wcor)
    total = a.weight + a.weight_correction
    newweight = w + wcor
    delta = (origin - a.origin) + (offset - a.offset)
    # Keep the heavier partial as anchor. A tiny-mass first outlier must
    # not erase small differences within the dominant population.
    crossweight = min(oldweight,newweight) * (max(oldweight,newweight) / total)
    a.norm2 += norm2 + delta * (crossweight * delta)
    if oldweight >= newweight
        a.offset += (newweight / total) * delta
    else
        a.origin = origin
        a.offset = offset - (oldweight / total) * delta
    end
end

function _merge_totals!(a, b)
    a.count += b.count
    a.volume, a.volume_correction = _compensated_add(a.volume, a.volume_correction, b.volume)
    a.volume, a.volume_correction = _compensated_add(a.volume, a.volume_correction, b.volume_correction)
    a.integral, a.correction = _compensated_add(a.integral, a.correction, b.integral)
    a.integral, a.correction = _compensated_add(a.integral, a.correction, b.correction)
    a.numerator, a.numerator_correction = _compensated_add(a.numerator, a.numerator_correction, b.numerator)
    a.numerator, a.numerator_correction = _compensated_add(a.numerator, a.numerator_correction, b.numerator_correction)
    a.minimum = min(a.minimum, b.minimum)
    a.maximum = max(a.maximum, b.maximum)
    _add_moments!(a, b.weight, b.origin, b.offset, b.norm2, b.weight_correction)
    a
end

# A lost nonzero product can give an exactly zero mean even for a constant
# positive field when both numerator and denominator carry a small scale.
@inline function _weighted_product(f, w)
    term = f * w
    iszero(term) && !iszero(f) && !iszero(w) &&
        throw(DomainError((f,w), "nonzero weighted product underflowed in Float64"))
    return term
end

function _accumulate!(a, b::SelectedBlock, values, densities)
    for k in axes(values,3), j in axes(values,2), i in axes(values,1)
        @inbounds f = Float64(values[i,j,k])
        isfinite(f) || throw(DomainError(f, "non-finite field value in block $(b.block)"))
        @inbounds dv = cell_volume(b, i, j, k)
        isfinite(dv) && dv > 0 || throw(DomainError(dv, "unrepresentable cell volume in block $(b.block)"))
        if isnothing(densities)
            w = dv
        else
            @inbounds rho = Float64(densities[i,j,k])
            isfinite(rho) && rho >= 0 || throw(DomainError(rho, "density must be finite and nonnegative in block $(b.block)"))
            w = rho * dv
            isfinite(w) && (rho == 0 || w > 0) || throw(DomainError(w, "unrepresentable cell mass in block $(b.block)"))
            a.numerator, a.numerator_correction = _compensated_add(a.numerator, a.numerator_correction, _weighted_product(f, w))
        end
        a.count += 1
        a.volume, a.volume_correction = _compensated_add(a.volume, a.volume_correction, dv)
        a.integral, a.correction = _compensated_add(a.integral, a.correction, _weighted_product(f, dv))
        a.minimum = min(a.minimum, f)
        a.maximum = max(a.maximum, f)
        _add_moments!(a, w, f, 0.0, 0.0)
    end
    a
end

# Only the caller executes this function. A compute task owns the returned
# arrays and never sees an HDF5 handle.
function _read_selection(s, b, f, rho)
    values = read_field(s, f, b.block; ranges=b.ranges)
    densities = if isnothing(rho)
        nothing
    elseif f == rho
        values
    else
        read_field(s, rho, b.block; ranges=b.ranges)
    end
    (b, values, densities)
end

function _accumulate_chunk(payload)
    total = CellTotals()
    for (b, values, densities) in payload
        _accumulate!(total, b, values, densities)
    end
    total
end

function _reduce_selections(mesh, blocks, f, rho, workers, chunk_blocks)
    s = mesh.snapshot
    total = CellTotals()
    if workers == 1
        for b in blocks
            block_total = _accumulate!(CellTotals(), _read_selection(s, b, f, rho)...)
            _merge_totals!(total, block_total)
        end
        return total
    end
    # Each wave owns at most workers × chunk_blocks selections. @sync joins
    # outstanding compute even if a read or a worker fails. Data-independent
    # merge order makes scheduling irrelevant at fixed partitioning.
    chunks = Iterators.partition(blocks, min(s.nblocks, chunk_blocks))
    for wave in Iterators.partition(chunks, workers)
        tasks = Task[]
        @sync for chunk in wave
            payload = [_read_selection(s, b, f, rho) for b in chunk]
            push!(tasks, Threads.@spawn _accumulate_chunk(payload))
        end
        for task in tasks
            _merge_totals!(total, fetch(task))
        end
    end
    total
end

"""
    summarize(mesh_or_snapshot, field; region=Region(), density=nothing, workers=1, chunk_blocks=1)

Summarize one cell-sampled field at the stored snapshot time. Values retain
their source frame; custom fields require the user's interpretation.
Returns `count`, `volume = Σ dV`, `integral = Σ f dV`, `minimum`, `maximum`,
`weight = Σ w`, `mean = Σ w f / Σ w`, and population
`variance = Σ w(f-mean)² / Σ w`. The integral always uses physical cell
volume. For volume weighting, w=dV. Passing a density selector explicitly
declares that field to be mass density in compatible code units and sets
w=ρ dV. No EOS or dimensional scaling is assumed. Results use Float64
accumulation. Mean and extrema have field units, variance has squared field
units, and the integral has field units times code length cubed.

Zero-density cells contribute to extrema, cell count and volume integral.
An empty selection, nonfinite field/density, negative density, zero total
weight or overflow in a Float64 numerator, weight, or central-moment sum
raises an exception, even when a rescaled calculation could return a finite
ratio. A nonzero field contribution that rounds to zero in a volume or mass
numerator also raises an exception. Subnormal arithmetic has reduced relative
precision. Face samples and
unsupported coverage are rejected. `provenance` retains source and selection,
field identity, weighting and compute configuration.

With multiple workers the caller reads all HDF5 hyperslabs serially. Compute
tasks each own at most `chunk_blocks` blocks of selected fields. The maximum
number of active tasks is bounded by Julia's available threads and `workers`.
Work is linear in selected cells plus a pass over block metadata; numerical
agreement across partitions is subject to floating-point rounding. Construct
`CartesianMesh(snapshot)` once to reuse verified coordinates for many calls.
"""
function summarize(mesh::CartesianMesh, field; region::Region=Region(), density=nothing,
                   workers::Int=1, chunk_blocks::Int=1)
    workers > 0 && chunk_blocks > 0 ||
        throw(ArgumentError("workers and chunk_blocks must be positive integers"))
    s = mesh.snapshot
    f = resolve_field(s, field)
    rho = isnothing(density) ? nothing : resolve_field(s, density)
    f.centering === :cell && (isnothing(rho) || rho.centering === :cell) ||
        throw(ArgumentError("reductions require cell-sampled fields; use raw access for other centerings"))
    nworkers = min(workers, Threads.nthreads(), s.nblocks)
    a = _reduce_selections(mesh, selectedblocks(mesh; region), f, rho, nworkers, chunk_blocks)
    a.count > 0 || throw(ArgumentError("the selected region contains no cell centers"))
    weight = a.weight + a.weight_correction
    volume = a.volume + a.volume_correction
    weight == 0 && throw(DomainError(weight, "the total selected weight must be positive"))
    integral = a.integral + a.correction
    numerator = isnothing(rho) ? integral : a.numerator + a.numerator_correction
    mean = numerator / weight
    variance = a.norm2 / weight
    all(isfinite, (volume, integral, weight, mean, variance)) ||
        throw(OverflowError("reduction exceeded the Float64 accumulation range"))
    provenance = (; path=s.path, time=s.time, cycle=s.cycle, field=f,
                    field_selector=field, density=rho, density_selector=density, region,
                    selection=:half_open_cell_center, geometry=:cartesian,
                    weighting=isnothing(rho) ? :volume : :mass,
                    workers=nworkers, requested_workers=workers, chunk_blocks)
    (; a.count, volume, integral, a.minimum, a.maximum, weight, mean, variance, provenance)
end

summarize(s::Snapshot, field; kwargs...) = summarize(CartesianMesh(s), field; kwargs...)
