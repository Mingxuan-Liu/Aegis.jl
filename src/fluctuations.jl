struct _GuideEnergy{C}
    reference::NTuple{3,Float64}
    direction::NTuple{3,Float64}
    scaled_direction::NTuple{3,Float64}
    direction_norm::Float64
end

function _guide_vector(v, name)
    (v isa Tuple || v isa AbstractVector) && length(v) == 3 && all(x -> x isa Real, v) ||
        throw(ArgumentError("$name needs three finite Cartesian components"))
    vec = ntuple(i -> Float64(v[i]), 3)
    all(isfinite, vec) || throw(ArgumentError("$name must be finite in Float64 code units"))
    return vec
end

# Compensate the product rounded before the first FMA. In particular an
# exactly parallel pair has no residual from independently rounded products.
@inline function _cross_term(a, b, c, d)
    cd = c*d
    fma(a, b, -cd) + fma(-c, d, cd)
end

@inline function _guide_projection(k::_GuideEnergy{:transverse}, x, y, z)
    a,b,c = k.scaled_direction
    (_cross_term(y,c,z,b), _cross_term(z,a,x,c), _cross_term(x,b,y,a))
end

@inline function _guide_projection(k::_GuideEnergy{:parallel}, x, y, z)
    a,b,c = k.scaled_direction
    xa, yb, zc = x*a, y*b, z*c
    # Recover small signed projections when the large Cartesian terms cancel.
    sum, corr = _compensated_add(xa, 0.0, yb)
    sum, corr = _compensated_add(sum, corr, zc)
    residual = fma(x,a,-xa) + fma(y,b,-yb) + fma(z,c,-zc)
    (sum + (corr + residual),)
end

@inline function (k::_GuideEnergy)(x, y, z)
    # Each difference is hi + lo. A small reference subtracted from a huge
    # field can disappear in hi, yet dominate its projected component.
    diffs = map((b,r) -> _compensated_add(b,0.0,-r), (x,y,z), k.reference)
    projection = _guide_projection(k, map(first,diffs)...)
    residuals = map(last, diffs)
    if any(!iszero, residuals)
        projection = map(+, projection, _guide_projection(k,residuals...))
    end
    amplitude = hypot(projection...) / k.direction_norm
    return _mhd_result((0.5 * amplitude)*amplitude, !iszero(amplitude))
end

@inline function (::_GuideEnergy)(energy)
    isfinite(energy) && energy >= 0 ||
        throw(DomainError(energy, "a stored magnetic fluctuation energy must be finite and nonnegative"))
    return energy
end

"""
    magnetic_fluctuation(snapshot, NewtonianMHD(), component; reference,
                         direction=reference, magnetic=nothing, stored=nothing)

Return a scalar field selector for parallel or transverse magnetic fluctuation
energy density in the Newtonian simulation frame (Cartesian axes, μ=1).
Subtract the uniform code-unit `reference` from the stored cell-centered B.
For unit direction n, `:parallel` gives ½(δB⋅n)² and `:transverse` gives
½|δB×n|², both in code energy per volume. The direction must be nonzero.

Supply three verified component selectors via `magnetic` (default Bcc1–3),
or declare one equivalent cell-centered scalar via `stored`. The latter
must use the given reference, direction, units and frame. Source and
declarations are retained in product metadata. See `docs/fluctuations.md`,
Stone et al. (2020) §3.1 Eq. (8), and Bell (2004) §2.
"""
function magnetic_fluctuation(s::Snapshot, ::NewtonianMHD, component::Symbol; reference,
                              direction=reference, magnetic=nothing, stored=nothing)
    s.coordinates == "cartesian" || throw(ArgumentError("a guide-field projection requires Cartesian components"))
    component in (:parallel, :transverse) || throw(ArgumentError("choose :parallel or :transverse"))
    r = _guide_vector(reference, "reference")
    d = _guide_vector(direction, "direction")
    largest = maximum(abs, d)
    iszero(largest) && throw(ArgumentError("supply a nonzero direction when using a zero reference"))
    # A common power of two preserves component ratios and keeps |d| finite.
    _, exponent = frexp(largest)
    g = map(x -> ldexp(x, -exponent), d)
    any(i -> !iszero(d[i]) && iszero(g[i]), 1:3) &&
        throw(ArgumentError("direction components exceed the supported Float64 dynamic range"))
    kernel = _GuideEnergy{component}(r, d, g, hypot(g...))
    name = "magnetic_fluctuation_$component"
    units = "code energy density"
    if !isnothing(stored)
        isnothing(magnetic) || throw(ArgumentError("choose stored energy or magnetic components"))
        _product_ref(s, stored) isa FieldRef || throw(ArgumentError("stored energy must select an on-disk cell field"))
        return cell_map(kernel, s, stored; name, units)
    end
    selectors = isnothing(magnetic) ? ("Bcc1","Bcc2","Bcc3") : Tuple(magnetic)
    length(selectors) == 3 || throw(ArgumentError("supply three Cartesian magnetic component selectors"))
    refs = map(f -> _product_ref(s, f), selectors)
    all(f -> f isa FieldRef, refs) && allunique(refs) ||
        throw(ArgumentError("magnetic components must select three distinct on-disk cell fields"))
    return cell_map(kernel, s, selectors...; name, units)
end

function _product_field_meta(m::CellMap{<:_GuideEnergy{C}}) where C
    (; kind=:magnetic_fluctuation, m.name, m.units, component=C,
       source=length(m.inputs) == 1 ? :stored : :reconstructed,
       reference=m.f.reference, direction=m.f.direction, inputs=m.sources,
       centering=:cell, frame=:simulation, fluid_model=:newtonian, permeability=1.0)
end
