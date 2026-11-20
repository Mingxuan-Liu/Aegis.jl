"""
    SpecialRelativisticMHD()

Declare ordinary flat-spacetime, Cartesian ideal SRMHD primitives from a
verified Athena++ producer. The stored velocity components are the spatial
four-velocity in a single inertial lab frame, in units with light speed one.
The cell-centered magnetic field is measured in the lab with magnetic
permeability one. This requires the run's fluid build and output convention.
See Stone et al. (2020), §4.1, and the pinned writer in `docs/relativity.md`.
"""
struct SpecialRelativisticMHD end

const _SR_LAB_VELOCITIES = (:lab_velocity_1, :lab_velocity_2, :lab_velocity_3)

@inline _sr_lorentz(x, y, z) = _mhd_result(hypot(1.0, x, y, z), true)

# Use u/W rather than 1/sqrt(1-v²), which loses the boost once v rounds to 1.
@inline function _sr_lab_component(i, x, y, z)
    u = (x, y, z)[i]
    return _mhd_result(u / _sr_lorentz(x, y, z), u != 0)
end

function _sr_magnetic_pressure(u1, u2, u3, b1, b2, b3)
    w = _sr_lorentz(u1, u2, u3)
    # Scale u and W by the same binary exponent. Independently rounding each
    # v_i before the dot product can fabricate a parallel field at high W.
    mantissa, power = frexp(w)
    x, y, z = ldexp(u1, -power), ldexp(u2, -power), ldexp(u3, -power)
    px, py, pz = x*b1, y*b2, z*b3
    sum12, err = _compensated_add(px, 0.0, py)
    sum123, err = _compensated_add(sum12, err, pz)
    prod_err = fma(x, b1, -px) + fma(y, b2, -py) + fma(z, b3, -pz)
    projection = (sum123 + (err + prod_err)) / mantissa
    bprime = hypot(b1/w, b2/w, b3/w, projection)
    return _mhd_result((0.5*bprime)*bprime, !iszero(b1) || !iszero(b2) || !iszero(b3))
end

function _sr_stored_value(name, value)
    valid = if name === :lorentz_factor
        value >= 1
    elseif name === :comoving_magnetic_pressure
        value >= 0
    else
        # A subluminal component can round to either endpoint on disk.
        abs(value) <= 1
    end
    valid || throw(DomainError(value, "stored $name is outside its physical range"))
    return value
end

# Locators may point to a custom producer's arbitrary names. In the pinned
# writer's native groups the standard roles are known and disjoint.
function _sr_input(s, field, name, role)
    ref = _mhd_input(s, field, name)
    standard_role = if ref.dataset in ("prim", "cons", "B", "hydro")
        if ref.name in ("vel1", "vel2", "vel3")
            :four_velocity
        elseif ref.name in ("Bcc1", "Bcc2", "Bcc3")
            :magnetic
        elseif ref.name in ("rho", "press", "dens", "mom1", "mom2", "mom3", "Etot")
            :other
        end
    end
    isnothing(standard_role) || standard_role === role ||
        throw(ArgumentError("stored $(ref.dataset)/$(ref.name) is incompatible with $role for $name"))
    return ref
end

function _sr_selectors(selectors, defaults, name, role)
    isnothing(selectors) && return defaults
    selectors isa Union{Tuple,AbstractVector} && length(selectors) == 3 ||
        throw(ArgumentError("supply three Cartesian $role field selectors for $name"))
    return Tuple(selectors)
end

"""
    quantity(snapshot, model::SpecialRelativisticMHD, name; four_velocity, magnetic, stored=nothing)

Select `:lorentz_factor` (dimensionless), `:lab_velocity_1` through
`:lab_velocity_3` (Cartesian components in units of c) or
`:comoving_magnetic_pressure` (fluid-frame code energy density). On each
cell the spatial four-velocity u gives W = √(1 + |u|²) and v = u/W. The
ideal-fluid magnetic invariant gives p′_B = (|B|²/W² + (v⋅B)²)/2.

Supply cell-centered spatial `four_velocity` (default vel1–vel3) and lab
`magnetic` field (default Bcc1–Bcc3) selectors in x1,x2,x3 order. By
setting `stored` to a single field selector the caller declares a saved
equivalent with exactly the requested units, centering and frame. Its values
are read directly and checked for physical range. Descriptors expose
`source` and `inputs` and work in `read_field`, `summarize` and cell
products. See Stone et al. (2020), §4.1, Mignone & Bodo (2006), Eqs. (4)–(6),
and `docs/relativity.md` for the output mapping and checked examples.
"""
function quantity(s::Snapshot, model::SpecialRelativisticMHD, name::Symbol;
                  four_velocity=nothing, magnetic=nothing, stored=nothing)
    s.coordinates == "cartesian" ||
        throw(ArgumentError("SpecialRelativisticMHD requires Cartesian inertial-frame output"))
    axis = findfirst(==(name), _SR_LAB_VELOCITIES)
    units = if name === :lorentz_factor
        "1"
    elseif name === :comoving_magnetic_pressure
        "code energy density"
    elseif !isnothing(axis)
        "c"
    else
        throw(ArgumentError("unsupported SR primitive diagnostic '$name'; see docs/relativity.md"))
    end
    if !isnothing(stored)
        isnothing(four_velocity) && isnothing(magnetic) ||
            throw(ArgumentError("choose either a stored equivalent or reconstruction input selectors"))
        inputs = (_sr_input(s, stored, name, :stored),)
        kernel = x -> _sr_stored_value(name, x)
        return MHDQuantity(s.path, String(name), units, :cell, model, inputs, kernel, :stored, (stored,))
    end
    us = _sr_selectors(four_velocity, ("vel1","vel2","vel3"), name, :four_velocity)
    u = map(f -> _sr_input(s, f, name, :four_velocity), us)
    if name === :comoving_magnetic_pressure
        bs = _sr_selectors(magnetic, ("Bcc1","Bcc2","Bcc3"), name, :magnetic)
        b = map(f -> _sr_input(s, f, name, :magnetic), bs)
        inputs = (u..., b...)
        selectors = (us..., bs...)
        kernel = _sr_magnetic_pressure
    else
        isnothing(magnetic) || throw(ArgumentError("$name uses only the spatial four-velocity"))
        inputs = u
        selectors = us
        kernel = isnothing(axis) ? _sr_lorentz : (x,y,z) -> _sr_lab_component(axis,x,y,z)
    end
    allunique(inputs) || throw(ArgumentError("supply distinct stored component fields for $name"))
    return MHDQuantity(s.path, String(name), units, :cell, model, inputs, kernel, :reconstructed, selectors)
end

function _product_field_meta(f::MHDQuantity{SpecialRelativisticMHD})
    frame = f.name == "comoving_magnetic_pressure" ? :fluid_rest : :inertial_lab
    return (; kind=:special_relativistic, f.source, f.name, f.units, frame, f.centering,
              velocity_convention=:spatial_four_velocity, light_speed=1.0,
              magnetic_permeability=1.0, inputs=map(_product_meta, f.selectors, f.inputs))
end
