module PeriodicSpectraExample

using Aegis
include("uniform_fixture.jl")

"""Write an analytic magnetic wave in a periodic cube of side π in code lengths."""
function write_wave(path)
    n = 8
    shape = (n, n, n)
    # Cell midpoints span [-π/2, π/2). Only x varies.
    wave = [sinpi(2*((i - 0.5)/n - 0.5)) for i in 1:n, j in 1:n, k in 1:n]
    density = (2 .+ wave).^2
    fields = (density, zeros(shape), zeros(shape), fill(2.0, shape),
              fill(3.0, shape), 2 .* wave, zeros(shape))
    names = ["rho", "vel1", "vel2", "vel3", "Bcc1", "Bcc2", "Bcc3"]
    write_uniform_snapshot(path, fields; lengths=(π, π, π), tiles=(2, 1, 1), names)
end

"""Analyze the manufactured wave in its simulation frame and in the gas rest frame."""
function analyze_wave(path)
    open_snapshot(path) do s
        mesh = CartesianMesh(s)
        periodic = (true, true, true)
        edges = 1.0:2.0:15.0  # angular wavenumbers, with centers 2, 4, …
        B = ("Bcc1", "Bcc2", "Bcc3")
        v = ("vel1", "vel2", "vel3")
        magnetic = shell_spectrum(mesh, B; periodic, edges)
        kinetic = shell_spectrum(mesh, v; density="rho", periodic, edges)
        centered = shell_spectrum(mesh, v; density="rho", periodic, edges, remove_mean=true)

        # The gas flows at 2 along x3 everywhere in this manufactured state.
        v3_rest = cell_map(v3 -> v3 - 2, s, "vel3"; name="v3 in gas rest frame")
        rest = shell_spectrum(mesh, ("vel1", "vel2", v3_rest); density="rho", periodic, edges)
        compensated = magnetic.k.^(5/3) .* magnetic.density
        (; magnetic, kinetic, centered, rest, compensated)
    end
end

function main()
    mktempdir() do dir
        r = analyze_wave(write_wave(joinpath(dir, "wave.athdf")))
        for (name, s) in (("magnetic", r.magnetic), ("kinetic", r.kinetic))
            closure = sum(s.shell_energy) + s.dc_energy + s.excluded_energy
            println(name, "  first shell=", round(s.shell_energy[1]; digits=8),
                    "  DC=", round(s.dc_energy; digits=8),
                    "  spectral sum=", round(closure; digits=8),
                    "  real-space energy=", round(s.real_energy; digits=8))
        end
        println("first compensated magnetic ordinate=", round(r.compensated[1]; digits=8))
        println("kinetic wave energy after centering FFT input=", round(sum(r.centered.shell_energy); digits=8))
        println("kinetic wave energy in gas rest frame=", round(sum(r.rest.shell_energy); digits=8))
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PeriodicSpectraExample.main()
end
