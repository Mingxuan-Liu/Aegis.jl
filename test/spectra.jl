include(joinpath(@__DIR__, "..", "examples", "uniform_fixture.jl"))

# Full direct DFT, one sum over all physical cells for every signed frequency.
# This oracle has no real-FFT plane weights, half grid, or FFT library.
function _direct_shells(fields, lengths, edges)
    dims = size(first(fields))
    n = prod(dims)
    frequencies = ntuple(d -> (-fld(dims[d],2)):(cld(dims[d],2)-1), 3)
    total = zeros(length(edges)-1)
    counts = zeros(Int, length(total))
    dc = excluded = 0.0
    outside_count = 0
    for v in Iterators.product(frequencies...)
        kmag = sqrt(sum((2π*v[d]/lengths[d])^2 for d in 1:3))
        power = 0.0
        for f in fields
            coeff = sum(ComplexF64(f[ix]) * cispi(-2sum((ix[d]-1)*v[d]/dims[d] for d in 1:3))
                        for ix in CartesianIndices(f)) / n
            power += abs2(coeff)/2
        end
        if all(iszero, v)
            dc = power
        elseif first(edges) <= kmag <= last(edges)
            bin = findfirst(i -> kmag < edges[i+1], eachindex(total))
            b = isnothing(bin) ? length(total) : bin
            total[b] += power
            counts[b] += 1
        else
            excluded += power
            outside_count += 1
        end
    end
    (; total, counts, dc, excluded, outside_count)
end

@testset "3D global Fourier shells and Parseval identity" begin
    periodic = (true,true,true)
    mktempdir() do dir
        for dims in ((5,4,3),(6,4,4))
            lengths = (10.0, 6.0, 3.5)
            # Includes nonzero DC, oblique modes, high spatial frequency,
            # and the x Nyquist plane when nx is even.
            a = [2.0 + 1.7cospi(2*(i/dims[1]+j/dims[2]-k/dims[3])) +
                 0.2cospi(2*fld(dims[1],2)*i/dims[1]) for
                 i in 0:dims[1]-1, j in 0:dims[2]-1, k in 0:dims[3]-1]
            b = [(-1)^j * cospi(2i/dims[1]) + 0.5 for
                 i in 0:dims[1]-1, j in 0:dims[2]-1, k in 0:dims[3]-1]
            rho = [1 + 0.3cospi(2j/dims[2]) for
                   i in 0:dims[1]-1, j in 0:dims[2]-1, k in 0:dims[3]-1]
            tiles = (iseven(dims[1]) ? 2 : 1, 2, iseven(dims[3]) ? 2 : 1)
            p = write_uniform_snapshot(joinpath(dir, "grid.athdf"), (a,b,rho); lengths, tiles)
            edges = [0.3, 1.3, 2.4, 3.1, 3.9] # incomplete shells deliberately
            open_snapshot(p) do s
                mesh = CartesianMesh(s)
                for weighted in (false,true), center in (false,true)
                    fields = weighted ? (sqrt.(rho).*a, sqrt.(rho).*b) : (a,b)
                    mu(f) = sum(f)/length(f)
                    input = center ? map(f -> f .- mu(f), fields) : fields
                    oracle = _direct_shells(input, lengths, edges)
                    spec = shell_spectrum(mesh, ("component1","component2");
                                          periodic, edges, centers=:geometric,
                                          density=weighted ? "component3" : nothing,
                                          remove_mean=center)
                    # Full DFT sums form complex phases numerically. A tolerance
                    # of 1e-13 also accommodates bins with only roundoff energy.
                    @test spec.mode_counts == oracle.counts
                    @test spec.excluded_modes == oracle.outside_count
                    @test sum(spec.mode_counts) + spec.excluded_modes + 1 == prod(dims)
                    @test spec.shell_energy ≈ oracle.total atol=1e-13
                    @test spec.density .* diff(edges) ≈ spec.shell_energy atol=1e-13
                    @test spec.dc_energy ≈ oracle.dc atol=1e-13
                    @test spec.excluded_energy ≈ oracle.excluded atol=1e-13
                    @test spec.real_energy ≈ sum(f -> sum(abs2,f), input)/(2prod(dims))
                    @test spec.real_energy ≈ sum(spec.shell_energy) + spec.excluded_energy + spec.dc_energy atol=1e-13
                    @test spec.original_mean_energy ≈ sum(f -> mu(f)^2, fields)/2
                    @test spec.k ≈ sqrt.(edges[1:end-1] .* edges[2:end])
                end
            end
        end
        dims = (8,3,3)
        lengths = (2π,2π,2π)
        wave = [4cospi(2*2i/8) for i in 0:7, j in 0:2, k in 0:2]
        nyquist = [3.0*(-1)^i for i in 0:7, j in 0:2, k in 0:2]
        guide = fill(5.0,dims)
        p = write_uniform_snapshot(joinpath(dir,"waves.athdf"), (wave,nyquist,guide); lengths,tiles=(2,1,3))
        open_snapshot(p) do s
            spec = shell_spectrum(s, ("component1","component2","component3");
                                  periodic, edges=[1.0,2.0,3.0,4.0])
            @test spec.shell_energy ≈ [0.,4.,4.5] atol=2e-14
            @test spec.dc_energy ≈ 12.5
            full = shell_spectrum(s, ("component1","component2","component3"); periodic)
            @test sum(full.shell_energy) ≈ 8.5
            @test full.excluded_energy == full.excluded_modes == 0
            @test full.real_energy ≈ 21
            incompatible = (interpret_field(s, "component2"; units="magnetic field"),
                            interpret_field(s, "component3"; units="velocity"))
            @test_throws ArgumentError shell_spectrum(s, incompatible; periodic)
            @test_throws ArgumentError shell_spectrum(s, (); periodic)
            @test_throws ArgumentError shell_spectrum(s, ("component1",); periodic=(true,true,false))
            @test_throws ArgumentError shell_spectrum(s, ("component1",); periodic, memory_budget=1)
            @test_throws ArgumentError shell_spectrum(s, ("component1",); periodic, edges=[0,1,1])
            @test_throws ArgumentError shell_spectrum(s, ("component1",); periodic, centers=:geometric)
            @test_throws DomainError shell_spectrum(s, ("component3",); periodic, density="component2")
        end
        h5open(p, "r+") do h
            centers = read(h["x1v"])
            centers[1,:] .+= 0.005 # still in cells, inconsistent sampling lattice
            h["x1v"][:,:] = centers
        end
        open_snapshot(p) do s
            @test_throws ArgumentError shell_spectrum(s, ("component1",); periodic)
        end
        open_snapshot(write_fixture(joinpath(dir, "nonuniform.athdf"))) do s
            @test_throws ArgumentError shell_spectrum(s, ("tracer",); periodic)
        end
        open_snapshot(write_fixture(joinpath(dir, "slab.athdf"); blocksize=(3,2,1))) do s
            @test_throws ArgumentError shell_spectrum(s, ("tracer",); periodic)
        end
    end
end
