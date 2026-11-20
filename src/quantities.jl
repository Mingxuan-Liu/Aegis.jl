"""Ideal adiabatic gas with finite, constant ratio of specific heats `gamma > 1`."""
struct IdealGas
    gamma::Float64
    function IdealGas(gamma::Real)
        g = Float64(gamma)
        isfinite(g) && g > 1 || throw(ArgumentError("ideal-gas gamma must be finite and > 1"))
        new(g)
    end
end

"""Isothermal gas with a positive, constant isothermal sound speed in fluid code units."""
struct Isothermal
    sound_speed::Float64
    function Isothermal(sound_speed::Real)
        c = Float64(sound_speed)
        isfinite(c) && c > 0 || throw(ArgumentError("isothermal sound speed must be finite and > 0"))
        new(c)
    end
end

"""
    NewtonianMHD(; eos=nothing, velocity_unit_cgs=nothing, mean_particle_mass_g=nothing)

Declare Newtonian gas primitives in the snapshot frame and Athena++ magnetic
units (permeability one). Choose `IdealGas(gamma)` or `Isothermal(c_iso)` for
EOS-dependent quantities. For kelvin temperatures, give the fluid velocity
unit in cm/s and the mean mass per thermal particle in grams. These scales
are independent of any artificial light speed used by particle dynamics.
"""
struct NewtonianMHD{E}
    eos::E
    velocity_unit_cgs::Union{Nothing,Float64}
    mean_particle_mass_g::Union{Nothing,Float64}
    function NewtonianMHD(; eos=nothing, velocity_unit_cgs=nothing, mean_particle_mass_g=nothing)
        eos isa Union{Nothing,IdealGas,Isothermal} ||
            throw(ArgumentError("choose IdealGas(gamma), Isothermal(c_iso), or omit the EOS"))
        v = isnothing(velocity_unit_cgs) ? nothing : Float64(velocity_unit_cgs)
        m = isnothing(mean_particle_mass_g) ? nothing : Float64(mean_particle_mass_g)
        for (name, value) in ((:velocity_unit_cgs,v), (:mean_particle_mass_g,m))
            isnothing(value) || (isfinite(value) && value > 0) ||
                throw(ArgumentError("$name must be positive and finite"))
        end
        new{typeof(eos)}(eos,v,m)
    end
end

struct MHDQuantity{M,I,F,S}
    path::String
    name::String
    units::String
    centering::Symbol
    model::M
    inputs::I
    kernel::F
    source::Symbol
    selectors::S
end

@inline function _gas_density(x)
    x > 0 || throw(DomainError(x, "Newtonian primitive gas density must be positive"))
    x
end

@inline function _gas_pressure(x)
    x >= 0 || throw(DomainError(x, "gas pressure must be nonnegative"))
    x
end

# All scalars reaching this guard are magnitudes, thermodynamic scalars or
# ratios. Preserve physical zeros while rejecting a false zero after underflow.
@inline function _mhd_result(value, nonzero::Bool)
    isfinite(value) && (!nonzero || value != 0) ||
        throw(DomainError(value, "local MHD diagnostic exceeds Float64 range; rescale the code units"))
    return value
end

@inline function _mhd_ratio(x, denominator)
    denominator > 0 || throw(DomainError(denominator, "zero denominator in local MHD ratio"))
    isfinite(denominator) && isfinite(x) ||
        throw(DomainError((x,denominator), "unrepresentable intermediate in local MHD ratio"))
    return _mhd_result(x/denominator, x != 0)
end

@inline _iso_pressure(r,c) = _mhd_result((_gas_density(r)*c)*c, true)

@inline _vecmag(x,y,z) = hypot(x,y,z)
@inline function _pmag(x,y,z)
    b = _vecmag(x,y,z)
    _mhd_result((0.5*b)*b, b != 0)
end

@inline _va(r,x,y,z) = _mhd_ratio(_vecmag(x,y,z), sqrt(_gas_density(r)))
@inline function _ekin(r,x,y,z)
    speed = _vecmag(x,y,z)
    v = sqrt(_gas_density(r))*speed
    return _mhd_result((0.5*v)*v, speed != 0)
end

function _mhd_recipe(name, model, rho, p, v, b)
    speed = "code velocity"
    energy = "code energy density"
    name === :speed && return (v, _vecmag, speed)
    name === :magnetic_magnitude && return (b, _vecmag, "code magnetic field")
    name === :kinetic_energy_density && return ((rho,v...), _ekin, energy)
    name in (:magnetic_pressure,:magnetic_energy_density) && return (b, _pmag, energy)
    name === :alfven_speed && return ((rho,b...), _va, speed)
    name === :alfven_mach && return ((rho,b...,v...),
        (r,b1,b2,b3,v1,v2,v3) -> _mhd_ratio(_vecmag(v1,v2,v3), _va(r,b1,b2,b3)), "1")

    eos = model.eos
    if name in (:gas_pressure,:plasma_beta)
        if eos isa Isothermal
            c = eos.sound_speed
            name === :gas_pressure && return ((rho,), r -> _iso_pressure(r,c), energy)
            return ((rho,b...), (r,x,y,z) -> _mhd_ratio(_iso_pressure(r,c), _pmag(x,y,z)), "1")
        end
        name === :gas_pressure && return ((p,), _gas_pressure, energy)
        return ((p,b...), (pg,x,y,z) -> _mhd_ratio(_gas_pressure(pg), _pmag(x,y,z)), "1")
    end

    name in (:thermal_energy_density,:sound_speed,:sonic_mach,:temperature) ||
        throw(ArgumentError("unknown Newtonian gas quantity '$name'; see docs/quantities.md"))
    isnothing(eos) && throw(ArgumentError("$name needs an explicit IdealGas(gamma) or Isothermal(c_iso) closure"))
    if name === :temperature
        (isnothing(model.mean_particle_mass_g) || isnothing(model.velocity_unit_cgs)) &&
            throw(ArgumentError("temperature in kelvin needs velocity_unit_cgs and mean_particle_mass_g"))
        # Exact Boltzmann constant in erg/K. P_unit/rho_unit = v_unit^2.
        scale = (model.mean_particle_mass_g/1.380649e-16) * model.velocity_unit_cgs^2
        isfinite(scale) && scale > 0 || throw(DomainError(scale, "unrepresentable kelvin conversion"))
        if eos isa Isothermal
            return ((), () -> _mhd_result((scale*eos.sound_speed)*eos.sound_speed,true), "K")
        end
        return ((p,rho), (pg,r) -> _mhd_result(_mhd_ratio(_gas_pressure(pg),_gas_density(r))*scale, pg != 0), "K")
    elseif eos isa Isothermal
        c = eos.sound_speed
        name === :sound_speed && return ((), () -> c, speed)
        name === :sonic_mach && return (v, (x,y,z) -> _mhd_ratio(_vecmag(x,y,z),c), "1")
        throw(ArgumentError("the isothermal closure does not evolve thermal energy density; an adiabatic IdealGas is required"))
    end
    gamma = eos.gamma
    name === :thermal_energy_density && return ((p,), pg -> _mhd_ratio(_gas_pressure(pg),gamma-1), energy)
    cs(pg,r) = _mhd_result(sqrt(gamma) * (sqrt(_gas_pressure(pg))/sqrt(_gas_density(r))),pg != 0)
    name === :sound_speed && return ((p,rho), cs, speed)
    return ((p,rho,v...), (pg,r,x,y,z) -> _mhd_ratio(_vecmag(x,y,z), cs(pg,r)), "1")
end

"""
    quantity(snapshot, model::NewtonianMHD, name::Symbol; density="rho", pressure="press",
             velocity=("vel1","vel2","vel3"), magnetic=("Bcc1","Bcc2","Bcc3"))

Describe one local Newtonian diagnostic on Cartesian cell-centered primitives.
Use this descriptor in `read_field`, `summarize`, and other cell analyses.
Only needed source fields are resolved and read. Selectors can be original
names, qualified dataset/index pairs, or stored field descriptors.

Supported names, units, equations and singularities are listed in
`docs/quantities.md`. For example, `:magnetic_energy_density` is |B|²/2 in
code energy per unit volume and `:alfven_speed` is |B|/√ρ in code velocity.
Thermodynamic quantities depend on the chosen EOS. Temperatures in kelvin
need explicit conversion scales. Ratios with zero denominator are rejected.
See Stone et al. (2020) §3.1 and Sun & Bai (2023) §2.2.
"""
function quantity(s::Snapshot, model::NewtonianMHD, name::Symbol; density="rho", pressure="press",
                  velocity=("vel1","vel2","vel3"), magnetic=("Bcc1","Bcc2","Bcc3"))
    s.coordinates == "cartesian" || throw(ArgumentError("Newtonian primitive quantities currently require Cartesian components"))
    v, b = Tuple(velocity), Tuple(magnetic)
    length(v) == 3 && length(b) == 3 || throw(ArgumentError("supply three Cartesian selectors for each vector"))
    selectors, kernel, units = _mhd_recipe(name, model, density, pressure, v, b)
    inputs = map(field -> _mhd_input(s, field, name), selectors)
    MHDQuantity(s.path, String(name), units, :cell, model, inputs, kernel, :reconstructed, selectors)
end

function quantity(::Snapshot, model, ::Symbol; kwargs...)
    throw(ArgumentError("unsupported gas model $(typeof(model)); declare NewtonianMHD or SpecialRelativisticMHD for the verified fluid primitives"))
end

function _mhd_input(s, field, name)
    ref = resolve_field(s, field)
    ref isa FieldRef && ref.centering === :cell ||
        throw(ArgumentError("$name requires stored cell-centered inputs; check selector $field"))
    return ref
end

const NewtonianQuantity = MHDQuantity{M} where {M<:NewtonianMHD}

function resolve_field(s::Snapshot, q::MHDQuantity)
    q.path == s.path || throw(ArgumentError("derived field belongs to a different snapshot path"))
    foreach(ref -> resolve_field(s,ref), q.inputs)
    q
end

@inline function _mhd_value(kernel::F, xs::Vararg{Any,N}) where {F,N}
    all(isfinite, xs) || throw(DomainError(xs, "nonfinite input to local MHD diagnostic"))
    y = kernel(map(Float64,xs)...)
    isfinite(y) || throw(DomainError(y, "local MHD diagnostic exceeded Float64 range"))
    y
end

# HDF5 hyperslabs acquire their native element types at runtime. Cross this
# function barrier once per block, so Julia specializes the inner loop for the
# kernel and concrete input-array tuple, including mixed Float32/64 inputs.
function _cell_slab(evaluate::E, kernel::F, arrays::Tuple) where {E,F}
    shape = axes(first(arrays))
    all(a -> axes(a) == shape, arrays) || throw(DimensionMismatch("cell fields must cover the same slab"))
    values = Array{Float64}(undef, size(first(arrays)))
    for i in eachindex(values)
        samples = map(a -> @inbounds(a[i]), arrays)
        @inbounds values[i] = evaluate(kernel, samples...)
    end
    return values
end

function read_field(s::Snapshot, q::MHDQuantity, block::Integer;
                    ranges=ntuple(d -> 1:s.blocksize[d], 3))
    _checkblock(s, block)
    resolve_field(s,q)
    if isempty(q.inputs)
        shape = length.(_slab_indices(s,ranges))
        return fill(_mhd_value(q.kernel),shape)
    end
    arrays = map(ref -> read_field(s,ref,block; ranges),q.inputs)
    return _cell_slab(_mhd_value, q.kernel, arrays)
end
