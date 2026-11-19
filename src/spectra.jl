import FFTW
using LinearAlgebra: mul!

"""
    spectrum_workspace(shape; blocksize, bins=max(1, maximum(shape) ÷ 2))
    spectrum_workspace(snapshot; bins=...)

Estimate bytes for a global Float64 3D real FFT, including one complex
half-grid, two Float64 input blocks, a one-real-grid allowance for FFTW
workspace, frequency axes and shell tables. Components reuse the same arrays
and plan. The estimate excludes Julia's runtime, the open mesh/HDF5 caches
and caller-defined field implementations. FFTW's internal workspace depends
on its plan, platform and thread count, so the estimate is a preflight check
and does not impose an operating-system memory limit.
"""
function spectrum_workspace(shape::NTuple{3,Int}; blocksize::NTuple{3,Int},
                            bins::Integer=max(1, maximum(shape) ÷ 2))
    all(>(0), shape) && all(>(0), blocksize) && bins > 0 ||
        throw(ArgumentError("grid, block and bin dimensions must be positive"))
    n = prod(big.(shape))
    half = (big(shape[1] ÷ 2) + 1) * shape[2] * shape[3]
    real_bytes, half_bytes = 8n, 16half
    plan_allowance = real_bytes
    block_bytes = 16prod(big.(blocksize))
    table_bytes = 16sum(big.(shape)) + 64big(bins) + 64
    estimated_bytes = real_bytes + half_bytes + plan_allowance + block_bytes + table_bytes
    return (; real_bytes, half_bytes, plan_allowance, block_bytes, table_bytes, estimated_bytes)
end

spectrum_workspace(s::Snapshot; kwargs...) =
    spectrum_workspace(s.rootsize; blocksize=s.blocksize, kwargs...)

function _uniform_lengths(mesh)
    s = mesh.snapshot
    all(>(1), s.rootsize) ||
        throw(ArgumentError("three-dimensional spectra require more than one cell on each axis"))
    return ntuple(3) do d
        lo, hi, _ = s.rootgrid[d]
        n = s.rootsize[d]
        length_d = hi - lo
        dx = length_d / n
        isfinite(dx) && dx > 0 || throw(ArgumentError("invalid length on FFT axis $d"))
        # Root attributes, faces and centers may be written at different precisions.
        coord_eps = max(eps(eltype(s.faces[d])), eps(eltype(s.centers[d])),
                        eps(float(eltype(s.attributes["RootGridX$d"]))))
        tol = min(4 * coord_eps * max(abs(lo), abs(hi), length_d), 0.002dx)
        for (segment_id, seg) in enumerate(mesh.axes[d])
            offset = (segment_id-1) * s.blocksize[d]
            for (i, face) in enumerate(seg.faces)
                abs(face - (lo + (offset+i-1)*dx)) <= tol ||
                    throw(ArgumentError("FFT requires uniformly spaced Cartesian faces on axis $d"))
            end
            for (i, center) in enumerate(seg.centers)
                abs(center - (lo + (offset+i-0.5)*dx)) <= tol ||
                    throw(ArgumentError("FFT requires midpoint samples on the uniform lattice on axis $d"))
            end
        end
        length_d
    end
end

# DFT mode magnitudes in radians per coordinate length. Only the first
# physical axis is compressed by FFTW.jl's rfft with dims=(1,2,3).
function _wave_axes(shape, lengths)
    ntuple(3) do d
        n = shape[d]
        stop = d == 1 ? n ÷ 2 : n-1
        step = 2π / lengths[d]
        [step * min(m, n-m) for m in 0:stop]
    end
end

function _shell_edges(edges, axes, nbins, centers)
    if isnothing(edges)
        kmax = sqrt(sum(abs2(maximum(axis)) for axis in axes))
        e = collect(range(0.0, nextfloat(kmax); length=nbins+1))
    else
        e = Float64.(collect(edges))
    end
    length(e) >= 2 && all(isfinite, e) && first(e) >= 0 &&
        all(x -> isfinite(x) && x > 0, diff(e)) ||
        throw(ArgumentError("supply strictly increasing finite nonnegative angular wavenumber edges"))
    if centers === :geometric
        first(e) > 0 || throw(ArgumentError("geometric shell centers require strictly positive edges"))
        k = [sqrt(e[i])*sqrt(e[i+1]) for i in 1:length(e)-1]
    elseif centers === :arithmetic
        k = e[1:end-1] .+ diff(e) ./ 2
    else
        throw(ArgumentError("shell centers must be :arithmetic or :geometric"))
    end
    return e, k
end

function _load_component!(grid, s, component, rho)
    for b in 1:s.nblocks
        destination = ntuple(3) do d
            offset = s.locations[d,b] * s.blocksize[d]
            (offset+1):(offset+s.blocksize[d])
        end
        block = read_field(s, component, b)
        target = view(grid, destination...)
        if isnothing(rho)
            copyto!(target, block)
        else
            density = component == rho ? block : read_field(s, rho, b)
            for i in eachindex(target, block, density)
                r = Float64(density[i])
                isfinite(r) && r >= 0 || throw(DomainError(r, "spectrum density must be finite and nonnegative in block $b"))
                target[i] = sqrt(r) * Float64(block[i])
            end
        end
    end
    all(isfinite, grid) || throw(DomainError(component, "nonfinite transformed field samples"))
    return grid
end

function _bin_modes!(energy, residual, counts, coeffs, axes, edges, n, nx, record_counts)
    kx, ky, kz = axes
    invn = inv(Float64(n))
    dc = 0.5 * abs2(coeffs[1] * invn)
    outside = outside_residual = 0.0
    excluded_count = 0
    @inbounds for l in eachindex(kz), j in eachindex(ky), i in eachindex(kx)
        (i == j == l == 1) && continue
        # The kx=0 and (even nx) kx=Nyquist planes still contain both signs
        # of ky,kz. Other kx planes stand for two full-grid modes.
        w = i == 1 || (iseven(nx) && i == length(kx)) ? 1 : 2
        kmag = sqrt(kx[i]^2 + ky[j]^2 + kz[l]^2)
        power = (0.5w) * abs2(coeffs[i,j,l] * invn)
        if first(edges) <= kmag <= last(edges)
            bin = min(searchsortedlast(edges, kmag), length(energy))
            energy[bin], residual[bin] = _compensated_add(energy[bin], residual[bin], power)
            record_counts && (counts[bin] += w)
        else
            outside, outside_residual = _compensated_add(outside, outside_residual, power)
            excluded_count += w
        end
    end
    return (; dc, outside=outside + outside_residual, excluded_count)
end

"""
    shell_spectrum(mesh_or_snapshot, components::Tuple; periodic, edges=nothing,
                   centers=:arithmetic, density=nothing, remove_mean=false,
                   memory_budget=2^30, fft_threads=1)

Compute the global, isotropic shell spectrum of one or more real, co-located
components on a verified uniform 3D Cartesian grid. Explicitly confirm periodic
boundary conditions in all directions with `periodic=(true,true,true)`.
For one scalar pass `(selector,)`. Magnetic components give magnetic energy
density under the Newtonian permeability-one convention. Passing a mass density
transforms `sqrt(rho) * component`; with three Newtonian velocities its norm
is the volume-averaged kinetic energy density.

Each returned `shell_energy` is the full-grid sum of one half the squared
modulus of FFT coefficients normalized by cell count. `density` is per unit
angular wavenumber and `k` contains representative shell centers. Bins use
left-closed intervals and include the last right edge. DC is always separated.
Default linear edges cover the Fourier-grid corner. Custom `edges` can set
an inscribed cutoff. `mode_counts` and `excluded_modes` count full-grid modes
once, independent of the number of components.

`real_energy` is the real-space half squared norm per cell after optional mean
subtraction. The closure identity is
`sum(shell_energy) + dc_energy + excluded_energy ≈ real_energy`.
`original_mean_energy` records the original zero-mode contribution, including
when subtracting the component means before the FFT. Returned provenance
includes field identities, grid and box lengths, selections and conventions.
See FFTW manual §4.8.6, Julia AbstractFFTs `rfft`, and `docs/spectra.md`.
"""
function shell_spectrum(mesh::CartesianMesh, components::Tuple; periodic,
                        edges=nothing, centers=:arithmetic, density=nothing,
                        remove_mean::Bool=false, memory_budget::Integer=2^30,
                        fft_threads::Int=1)
    periodic === (true,true,true) ||
        throw(ArgumentError("global spectra require confirmed periodic boundaries on all three axes"))
    isempty(components) && throw(ArgumentError("supply at least one component selector"))
    fft_threads >= 1 || throw(ArgumentError("fft_threads must be positive"))
    s = mesh.snapshot
    refs = map(f -> resolve_field(s, f), components)
    declared_units = unique([f.units for f in components
                             if hasproperty(f, :units) && !isnothing(f.units)])
    length(declared_units) <= 1 ||
        throw(ArgumentError("component unit labels conflict; convert the fields to a common unit first"))
    rho = isnothing(density) ? nothing : resolve_field(s, density)
    all(f -> f.centering === :cell, refs) && (isnothing(rho) || rho.centering === :cell) ||
        throw(ArgumentError("spectrum components and density must be sampled at cell centers"))
    lengths = _uniform_lengths(mesh)
    shape = s.rootsize
    nbins = isnothing(edges) ? max(1, maximum(shape) ÷ 2) : length(edges)-1
    budget = spectrum_workspace(s; bins=nbins)
    budget.estimated_bytes <= memory_budget || throw(ArgumentError(
        "global FFT estimate $(budget.estimated_bytes) bytes exceeds memory_budget=$memory_budget"))
    waves = _wave_axes(shape, lengths)
    e, k = _shell_edges(edges, waves, nbins, centers)
    grid = Array{Float64}(undef, shape)
    half_shape = (shape[1] ÷ 2 + 1, shape[2], shape[3])
    coeffs = Array{ComplexF64}(undef, half_shape)
    plan = FFTW.plan_rfft(grid, (1,2,3); flags=FFTW.ESTIMATE, num_threads=fft_threads)
    energy = zeros(nbins)
    residual = zeros(nbins)
    counts = zeros(Int, nbins)
    n = length(grid)
    dc = excluded = real_energy = mean_energy = 0.0
    excluded_modes = 0
    for (q, field) in enumerate(refs)
        _load_component!(grid, s, field, rho)
        mu = sum(grid) / n
        mean_energy += 0.5mu^2
        remove_mean && (grid .-= mu)
        real_energy += 0.5*(sum(abs2, grid)/n)
        mul!(coeffs, plan, grid)
        part = _bin_modes!(energy, residual, counts, coeffs, waves, e, n, shape[1], q == 1)
        dc += part.dc
        excluded += part.outside
        excluded_modes = part.excluded_count
    end
    energy .+= residual
    spectral_density = energy ./ diff(e)
    all(isfinite, (real_energy, dc, excluded, mean_energy)) &&
        all(isfinite, spectral_density) || throw(OverflowError("spectrum exceeded Float64 dynamic range"))
    provenance = (; path=s.path, time=s.time, cycle=s.cycle, shape, lengths, components=refs,
                    selectors=components, density=rho, density_selector=density, remove_mean,
                    periodic, centers, wavenumber=:angular, fft_threads, workspace=budget)
    return (; edges=e, k, shell_energy=energy, density=spectral_density, mode_counts=counts,
              dc_energy=dc, excluded_energy=excluded, excluded_modes, real_energy,
              original_mean_energy=mean_energy, provenance)
end

shell_spectrum(s::Snapshot, components::Tuple; kwargs...) =
    shell_spectrum(CartesianMesh(s), components; kwargs...)
