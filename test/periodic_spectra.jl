include(joinpath(@__DIR__, "..", "examples", "periodic_spectra.jl"))

@testset "Executable periodic wave example" begin
    mktempdir() do dir
        r = PeriodicSpectraExample.analyze_wave(
            PeriodicSpectraExample.write_wave(joinpath(dir, "wave.athdf")))
        # A sinusoid of amplitude 2 has energy 1. Computation here needs only
        # sinpi, 8-point FFTs and Float64 reductions, so allow tens of ulps.
        for (s, dc) in ((r.magnetic, 4.5), (r.kinetic, 8.0), (r.centered, 0.0))
            @test s.k[1] == 2
            @test s.shell_energy[1] ≈ 1 atol=5e-14
            @test s.density[1] ≈ 0.5 atol=5e-14
            @test sum(s.shell_energy[2:end]) < 5e-28
            @test s.dc_energy ≈ dc atol=5e-14
            @test s.excluded_modes == 0
            @test s.real_energy ≈ 1 + dc atol=5e-14
            @test sum(s.shell_energy) + s.dc_energy + s.excluded_energy ≈ s.real_energy atol=5e-14
        end
        @test first(r.compensated) ≈ 2.0^(2/3) atol=5e-14
        @test r.centered.original_mean_energy ≈ 8 atol=5e-14
        @test all(iszero, r.rest.shell_energy)
        @test r.rest.real_energy == 0
    end
end
