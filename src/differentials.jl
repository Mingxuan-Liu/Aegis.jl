"""
    DifferentialGrid(mesh::CartesianMesh; boundary)

Prepare a complete level-zero Cartesian mesh for reconstructed first derivatives.
Supply three boundary modes. Use `:periodic` for a periodic pair, `:onesided`
for an interior quadratic fit at physical boundaries, and `:invariant` for a
physically homogeneous singleton direction. Other axes require at least three
centers. Reuse the grid inside the original snapshot's open lifetime.
"""
struct DifferentialGrid{M}
    mesh::M
    boundary::NTuple{3,Symbol}
    centers::NTuple{3,Vector{Float64}}
    blocks::Array{Int,3}
    function DifferentialGrid(mesh::CartesianMesh; boundary=nothing)
        boundary isa NTuple{3,Symbol} || throw(ArgumentError(
            "supply three axis boundary modes: :periodic, :onesided or :invariant"))
        s = mesh.snapshot
        for (n, mode) in zip(s.rootsize, boundary)
            if n == 1
                mode == :invariant || throw(ArgumentError(
                    "a singleton axis needs an explicit :invariant assumption"))
            elseif n < 3 || !(mode in (:periodic, :onesided))
                throw(ArgumentError("an active derivative axis needs ≥3 centers and :periodic or :onesided"))
            end
        end
        blocks = Array{Int}(undef, s.rootsize .÷ s.blocksize)
        for b in 1:s.nblocks
            loc = ntuple(d -> s.locations[d,b] + 1, 3)
            blocks[loc...] = b
        end
        centers = ntuple(d -> reduce(vcat, (a.centers for a in mesh.axes[d])), 3)
        new{typeof(mesh)}(mesh, boundary, centers, blocks)
    end
end

# Each term describes one directed partial derivative. A curl component uses
# a pair of signed terms; divergence uses three. Selectors retain unit labels.
struct DifferentialField{G,T}
    grid::G
    terms::T
    name::String
    units::Union{String,Nothing}
    centering::Symbol
    path::String
end

function _derivative_field(grid, specs)
    s = grid.mesh.snapshot
    terms = map(specs) do (field, axis, sign)
        f = resolve_field(s, field)
        f.centering === :cell || throw(ArgumentError("derivatives require cell-sampled fields"))
        units = hasproperty(field, :units) ? field.units : nothing
        (; field, source=f, axis, sign, units)
    end
    declared = unique([t.units for t in terms if !isnothing(t.units)])
    length(declared) <= 1 || throw(ArgumentError("derivative terms have incompatible field unit labels"))
    units = all(t -> !isnothing(t.units), terms) ? "($(only(declared)))/code_length" : nothing
    labels = [string(t.sign < 0 ? "-" : (i == 1 ? "" : "+"), "d", t.source.name, "/dx", t.axis)
              for (i,t) in enumerate(terms)]
    DifferentialField(grid, terms, join(labels), units, :cell, s.path)
end

"""
    gradient(grid, field)

Return three selectable first derivatives of a scalar, ordered by (x1,x2,x3).
Inputs and outputs are cell-centered, with output units of field per code length.
The three-center quadratic stencil uses neighboring native blocks and the grid's
declared physical boundary rule. See Fornberg (1988), doi:10.1090/S0025-5718-1988-0935077-0.
"""
gradient(grid::DifferentialGrid, field) = ntuple(d -> _derivative_field(grid, ((field,d,1),)), 3)

"""
    vector_gradient(grid, fields)

Return three rows of derivative selectors for the ordered Cartesian vector
components. Entry `[i][j]` is ∂fields[i]/∂xⱼ in component units per code length.
For cell-centered magnetic components this reconstructs their gradient tensor.
"""
vector_gradient(grid::DifferentialGrid, fields::Tuple{Any,Any,Any}) =
    map(f -> gradient(grid, f), fields)

"""
    divergence(grid, fields)

Return a selectable reconstructed cell-centered ∑ᵢ ∂fields[i]/∂xᵢ for a
Cartesian vector tuple in compatible units. Its units are vector units per
code length. The native constrained-transport divergence would require face
fluxes. See Stone et al. (2020), §3.2, doi:10.3847/1538-4365/ab929b.
"""
divergence(grid::DifferentialGrid, fields::Tuple{Any,Any,Any}) =
    _derivative_field(grid, ntuple(d -> (fields[d], d, 1), 3))

"""
    curl(grid, fields)

Return the three cell-centered components of ∇×fields in right-handed Cartesian
order as selectors. For velocity in code units these are the components of
vorticity in inverse code time. Stencils and boundary choices follow `gradient`.
"""
function curl(grid::DifferentialGrid, fields::Tuple{Any,Any,Any})
    ntuple(3) do i
        j, k = mod1(i+1,3), mod1(i+2,3)
        _derivative_field(grid, ((fields[k],j,1), (fields[j],k,-1)))
    end
end

function resolve_field(s::Snapshot, f::DifferentialField)
    s === f.grid.mesh.snapshot || throw(ArgumentError("derivative grid belongs to another snapshot handle"))
    return f
end

function _product_field_meta(f::DifferentialField)
    ingredients = map(f.terms) do t
        (; source=_product_meta(t.field,t.source), axis=t.axis, sign=t.sign)
    end
    (; kind=:reconstructed_derivative, name=f.name, units=something(f.units,"unspecified"), path=f.path,
       centering=f.centering, boundary=f.grid.boundary, method="three_center_quadratic", ingredients)
end

# Gather a single-axis strip into stencil order. A virtual index of zero or
# N+1 denotes one periodic halo plane. Fetch maximal contiguous ranges from
# each contributing MeshBlock, including blocks narrower than the stencil.
function _derivative_strip(grid, source, b, ranges, axis, firstrow, lastrow)
    s = grid.mesh.snapshot
    n = s.rootsize[axis]
    width = s.blocksize[axis]
    target_shape = ntuple(d -> d == axis ? lastrow-firstrow+1 : length(ranges[d]), 3)
    strip = Array{Float64}(undef, target_shape)
    row = firstrow
    while row <= lastrow
        globalcell = mod1(row,n)
        pos, offset = divrem(globalcell-1,width)
        run = min(width-offset, lastrow-row+1)
        loc = ntuple(d -> d == axis ? pos+1 : s.locations[d,b]+1, 3)
        slab_ranges = Base.setindex(ranges, (offset+1):(offset+run), axis)
        slab = read_field(s, source, grid.blocks[loc...]; ranges=slab_ranges)
        all(isfinite, slab) || throw(DomainError(source.name, "non-finite derivative stencil sample"))
        target = (row-firstrow+1):(row-firstrow+run)
        selectdim(strip, axis, target) .= slab
        row += run
    end
    return strip
end

function _center_gap(grid, axis, lower)
    xs = grid.centers[axis]
    if 1 <= lower < length(xs)
        return xs[lower+1] - xs[lower]
    end
    # Avoid loss of precision from adding the full period to the first center.
    segments = grid.mesh.axes[axis]
    return (last(last(segments).faces) - last(xs)) + (first(xs) - first(first(segments).faces))
end

function _add_partial!(out, grid, term, block, ranges)
    axis = term.axis
    bc = grid.boundary[axis]
    bc === :invariant && return
    s = grid.mesh.snapshot
    n = s.rootsize[axis]
    offset = s.blocksize[axis] * s.locations[axis,block]
    lo, hi = first(ranges[axis]) + offset, last(ranges[axis]) + offset
    left = bc === :periodic ? lo-1 : clamp(lo-1, 1, n-2)
    right = bc === :periodic ? hi+1 : clamp(hi+1, 3, n)
    strip = _derivative_strip(grid, term.source, block, ranges, axis, left, right)
    for (p, globalcell) in enumerate(lo:hi)
        a = bc === :periodic ? globalcell-1 : clamp(globalcell-1, 1, n-2)
        h1, h2 = _center_gap(grid, axis, a), _center_gap(grid, axis, a+1)
        # Quadratic derivative via the two consecutive secant slopes. Using
        # sample differences makes the derivative of a constant exactly zero.
        x_minus_b = globalcell == a ? -h1 : (globalcell == a+1 ? 0.0 : h2)
        blend = (2*x_minus_b + h1) / (h1+h2)
        ia = a - left + 1
        fa, fb, fc = ntuple(i -> selectdim(strip, axis, ia+i-1), 3)
        dest = selectdim(out, axis, p)
        for idx in eachindex(dest,fa,fb,fc)
            s1 = (fb[idx]-fa[idx])/h1
            s2 = (fc[idx]-fb[idx])/h2
            dest[idx] += term.sign * (s1 + blend*(s2-s1))
        end
    end
end

function read_field(s::Snapshot, field::DifferentialField, block::Integer;
                    ranges=ntuple(d -> 1:s.blocksize[d], 3))
    resolve_field(s, field)
    _checkblock(s, block)
    indices = _slab_indices(s,ranges)
    out = zeros(Float64, length.(indices))
    for term in field.terms
        _add_partial!(out, field.grid, term, Int(block), indices)
    end
    all(isfinite, out) || throw(OverflowError("non-finite reconstructed derivative"))
    return out
end
