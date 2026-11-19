# Measure the last native state of each supplied, single-stream resolution run.
# The output directory receives compact per-run HDF5 diagnostics for compare.py.
using Aegis
using HDF5
using SHA

function measure_run(directory, destination)
    deck_path = only(filter(p -> startswith(basename(p), "athinput"), readdir(directory; join=true)))
    deck = read_athinput(deck_path)
    periodic = ntuple(d -> all(side -> parameter(deck, "mesh", "$(side)x$(d)_bc") == "periodic",
                               ("i", "o")), 3)
    id = parameter(deck, "job", "problem_id")
    eta = parameter(deck, "problem", "eta_ohm", Float64)
    nu = parameter(deck, "problem", "nu_iso", Float64)
    records = only(values(discover_snapshots(directory)))
    path = only(r.path for r in records if r.time == last(records).time)
    open_snapshot(path) do snapshot
        shape = snapshot.rootsize
        lengths = ntuple(d -> snapshot.rootgrid[d][2]-snapshot.rootgrid[d][1], 3)
        all(==(lengths[1]), lengths) && all(==(shape[1]), shape) ||
            error("This publication comparison requires a cubic box and cubic grid")
        # Figure recipe from turbulence_shell_edges in the producing project.
        k0 = 2pi / lengths[1]
        kmax = sqrt(sum((k0 * (n ÷ 2))^2 for n in shape)) * 1.00001
        nbins = ceil(Int, 4log10(kmax/k0))
        edges = 10.0 .^ range(log10(k0), log10(kmax); length=nbins+1)
        mesh = CartesianMesh(snapshot)
        budget = 6 * 2^30
        result = @timed shell_spectrum(mesh, ("Bcc1","Bcc2","Bcc3");
                                      periodic, edges, centers=:geometric, memory_budget=budget)
        spectrum = result.value
        @assert isapprox(sum(spectrum.shell_energy) + spectrum.dc_energy + spectrum.excluded_energy,
                         spectrum.real_energy; rtol=1e-10)
        @assert sum(spectrum.mode_counts) == prod(shape)-1
        @assert iszero(spectrum.excluded_energy)
        out = joinpath(destination, splitext(basename(path))[1] * ".h5")
        h5open(out, "w") do h
            a = attributes(h)
            for (key, value) in pairs((; run=id, source=basename(path), t_code=snapshot.time,
                                     cycle=snapshot.cycle, n=first(shape), L=first(lengths),
                                     source_bytes=filesize(path), eta, nu,
                                     deck_sha256=bytes2hex(sha256(read(deck_path))),
                                     dc_energy=spectrum.dc_energy, real_energy=spectrum.real_energy,
                                     excluded_energy=spectrum.excluded_energy))
                a[string(key)] = value
            end
            for key in (:edges, :k, :shell_energy, :density, :mode_counts)
                h[string(key)] = getproperty(spectrum, key)
            end
        end
        println(basename(path), "  shells=", nbins, "  spectrum_seconds=", result.time,
                "  process_peak_bytes=", Int(Sys.maxrss()), "  real_energy=", spectrum.real_energy)
    end
end

length(ARGS) >= 2 || error("Usage: julia --project=. measure.jl DESTINATION RUN_DIRECTORY ...")
mkpath(first(ARGS))
foreach(dir -> measure_run(dir, first(ARGS)), ARGS[2:end])
