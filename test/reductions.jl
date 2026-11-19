@testset "Cartesian cell measures and weighted statistics" begin
    mktempdir() do dir
        for order in ([1,2], [2,1]), nz in (1, 2)
            p = write_fixture(joinpath(dir, "cells.athdf"); order, blocksize=(3,2,nz))
            open_snapshot(p) do s
                mesh = CartesianMesh(s)
                # x1 widths 1,...,6. The transverse measure is 3 or 9.
                area = nz == 1 ? 3 : 9
                for workers in (1, 2, 4), chunk_blocks in (1, 2, 3)
                    v = summarize(mesh, "tracer"; workers, chunk_blocks)
                    @test v.count == 6*2*nz
                    @test v.volume == v.weight == 21area
                    @test v.integral == 91area
                    @test (v.minimum, v.maximum) == (1, 6)
                    @test v.mean ≈ 91/21 rtol=2e-14
                    # E[x²]=Σi³/Σi; subtract (Σi²/Σi)² analytically.
                    @test v.variance ≈ 20/9 rtol=2e-14
                    m = summarize(mesh, "tracer"; density="rho", workers, chunk_blocks)
                    @test m.weight == 112area
                    @test m.integral == v.integral
                    @test m.mean ≈ 4.75 rtol=2e-14
                    @test m.variance ≈ 1.6875 rtol=2e-14
                    @test m.provenance.weighting == :mass
                    @test m.provenance.field.name == "tracer"
                    @test m.provenance.density.name == "rho"
                end
                # Centers at 2 and 4.5 are included, the center at 8 is not.
                # A tiny interval around the center at 2 selects its full cell.
                box = Region(x1=(2,8))
                sub = summarize(mesh, "tracer"; region=box)
                @test sub.volume == 5area
                @test sub.mean ≈ 13/5
                @test sub.variance ≈ 6/25
                block = only(selectedblocks(mesh; region=box))
                @test block.block == findfirst(==(1), order)
                @test block.ranges == (2:3, 1:2, 1:nz)
                @test cell_volume(block, 1, 2, nz) == 2*2*nz
                @test summarize(mesh, "tracer"; region=Region(x1=(2,2.01))).volume == 2area
                @test isempty(collect(selectedblocks(mesh; region=Region(x2=(10,11)))))
                @test_throws ArgumentError summarize(mesh, "tracer"; region=Region(x1=(2,2)))
                @test_throws ArgumentError summarize(mesh, "rho"; workers=0)
                @test_throws ArgumentError summarize(mesh, "rho"; chunk_blocks=0)
                @test_throws ArgumentError summarize(mesh, "rho"; density="absent")
            end
        end
        @test_throws ArgumentError Region(x2=(1,0))
        @test_throws ArgumentError Region(x3=(NaN,1))

        # More blocks than workers and a final short chunk exercise wave
        # boundaries. Reference uses exact integer weights and rational
        # centered deviations, independent of the online production algorithm.
        n = 33
        xbar = sum(BigInt(i)^2 for i in 1:n) // sum(1:n)
        exact_var = Float64(sum(i*(i-xbar)^2 for i in 1:n) / sum(1:n))
        p = write_fixture(joinpath(dir,"waves.athdf"); order=collect(11:-1:1), fields=("tracer",))
        h5open(p,"r+") do h
            a = read(h["extras"])
            HDF5.delete_object(h, "extras")
            h["extras"] = Float64.(a) .+ 2.0^48
        end
        open_snapshot(p) do s
            mesh = CartesianMesh(s)
            for (w, c) in ((1,1), (2,1), (3,2), (3,4))
                shifted = summarize(mesh, "tracer"; workers=w, chunk_blocks=c)
                @test shifted.mean == 2.0^48 + Float64(xbar)
                @test shifted.variance ≈ exact_var rtol=8e-14
            end
        end
    end
end

@testset "reject unsupported coverage and ambiguous centering" begin
    function invalid_mesh(change)
        mktempdir() do dir
            p = write_fixture(joinpath(dir,"bad.athdf"))
            h5open(change, p, "r+")
            open_snapshot(s -> CartesianMesh(s), p)
        end
    end
    function replace_attr(h, key, v)
        HDF5.delete_attribute(h, key)
        attributes(h)[key] = v
    end
    @test_throws ArgumentError invalid_mesh(h -> replace_attr(h, "Coordinates", "cylindrical"))
    @test_throws ArgumentError invalid_mesh(h -> replace_attr(h, "MaxLevel", Int32(1)))
    # A sliced/summed or ghosted block can have a perfectly valid data shape.
    @test_throws ArgumentError invalid_mesh(h -> replace_attr(h, "RootGridSize", Int32[10,2,2]))
    @test_throws ArgumentError invalid_mesh(h -> replace_attr(h, "RootGridSize", Int32[9,2,2]))
    @test_throws ArgumentError invalid_mesh(h -> h["LogicalLocations"][1,1] = 0)
    @test_throws ArgumentError invalid_mesh(h -> h["LogicalLocations"][1,1] = 2)
    @test_throws ArgumentError invalid_mesh(h -> h["x1f"][1,1] = 6.25) # interface gap
    @test_throws ArgumentError invalid_mesh(h -> h["x1f"][4,1] = 23.)  # root endpoint
    @test_throws ArgumentError invalid_mesh(h -> h["x2f"][2,1] = 1.1) # transverse mismatch
    @test_throws ArgumentError invalid_mesh(h -> h["x1v"][1,1] = 13.) # center outside
    @test_throws ArgumentError invalid_mesh(h -> h["x1f"][2,1] = 5.)  # negative width
    mktempdir() do dir
        p = write_fixture(joinpath(dir, "face.athdf"); fields=("B1","rho"))
        open_snapshot(p) do s
            @test_throws ArgumentError summarize(s,"B1")
            @test_throws ArgumentError summarize(s,"rho"; density="B1")
        end
    end
end

@testset "invalid values and zero density" begin
    mktempdir() do dir
        p = write_fixture(joinpath(dir,"weights.athdf"))
        function change_and_reduce(change; weighted=true)
            h5open(change, p, "r+")
            open_snapshot(p) do s
                summarize(s,"tracer"; density=weighted ? "rho" : nothing)
            end
        end
        # Set all mass in one cell only. Zero-weight cells are still included
        # in the volume integral, extrema, and count.
        m = change_and_reduce() do h
            gas = zeros(Float32,size(h["gas"]))
            gas[2,1,1,2,1] = 3
            h["gas"][:,:,:,:,:] = gas
        end
        @test m.count == 24
        @test (m.minimum, m.maximum, m.integral) == (1,6,819)
        @test (m.weight, m.mean, m.variance) == (6,2,0)
        @test_throws DomainError change_and_reduce(h -> h["gas"][2,1,1,2,1] = 0)
        @test_throws DomainError change_and_reduce(h -> h["gas"][2,1,1,2,1] = -1)
        @test_throws DomainError change_and_reduce(h -> h["gas"][2,1,1,2,1] = Inf32)
        @test_throws DomainError change_and_reduce(h -> h["extras"][1,1,1,1,2] = NaN32; weighted=false)
        # Explicit constant and zero fields.
        for constant in (3.5f0, 0f0)
            v = change_and_reduce(; weighted=false) do h
                h["extras"][:,:,:,:,2] = fill(constant,3,2,2,2)
            end
            @test (v.integral, v.mean, v.variance) == (189constant, constant, 0)
        end
    end
end

@testset "signed cancellation and tiny-mass outliers" begin
    mktempdir() do dir
        for order in ([1,2], [2,1])
            p = write_fixture(joinpath(dir,"precision.athdf"); blocksize=(2,1,1),
                              fields=("rho","tracer"), order)
            fvals = [1e16,1.,1.,-1e16]
            for extreme_weight in (false,true)
                values = extreme_weight ? [1e20,1.,2.,3.] : fvals
                weights = extreme_weight ? [1e-50,1.,1.,1.] : ones(4)
                h5open(p,"r+") do h
                    HDF5.delete_attribute(h,"RootGridX1")
                    attributes(h)["RootGridX1"] = [0.,4.,1.]
                    field = zeros(Float64,2,1,1,2,1)
                    rho = similar(field)
                    for (b,position) in pairs(order)
                        start = 2 * (position-1)
                        h["x1f"][:,b] = Float64[start,start+1,start+2]
                        h["x1v"][:,b] = [start+0.5,start+1.5]
                        field[:,1,1,b,1] = values[start+1:start+2]
                        rho[:,1,1,b,1] = weights[start+1:start+2]
                    end
                    HDF5.delete_object(h,"extras")
                    HDF5.delete_object(h,"gas")
                    h["extras"], h["gas"] = field, rho
                end
                w = BigFloat.(weights)
                x = BigFloat.(values)
                μ = sum(w .* x)/sum(w)
                v = sum(w .* (x .- μ).^2)/sum(w)
                open_snapshot(p) do s
                    mesh = CartesianMesh(s)
                    for (workers,nchunk) in ((1,1),(2,1),(2,typemax(Int)))
                        m = summarize(mesh,"tracer"; workers, chunk_blocks=nchunk, density="rho")
                        @test m.mean ≈ Float64(μ) rtol=2e-14
                        @test m.variance ≈ Float64(v) rtol=2e-14
                        if !extreme_weight
                            plain = summarize(mesh,"tracer"; workers, chunk_blocks=nchunk)
                            @test plain.mean == plain.integral/plain.weight == 0.5
                        end
                    end
                end
            end
        end
    end
end

@testset "small positive cell masses accumulate" begin
    mktempdir() do dir
        p = write_fixture(joinpath(dir, "light.athdf"); fields=("rho","tracer"),
                          order=[1], blocksize=(32768,1,1))
        n = 32768
        h5open(p,"r+") do h
            h["extras"][:,:,:,:,:] = ones(Float32, n,1,1,1,1)
            HDF5.delete_object(h, "gas")
            # Fixture width is i. Choose densities that give mass 1 to
            # the first cell and mass ~1e-16 to each remaining cell.
            dens = [i == 1 ? 1.0 : 1e-16 / i for i in 1:n]
            h["gas"] = reshape(dens, n,1,1,1,1)
        end
        open_snapshot(p) do s
            r = summarize(s, "tracer"; density="rho")
            @test r.weight ≈ 1.0 + (n-1)*1e-16 rtol=eps(Float64)
            @test r.mean == r.minimum == r.maximum == 1.0
            @test r.variance == 0.0
        end
    end
end

@testset "unrepresentable weighted numerator" begin
    mktempdir() do dir
        p = write_fixture(joinpath(dir,"tiny.athdf"); fields=("rho","tracer"),
                          order=[1], blocksize=(1,1,1))
        h5open(p,"r+") do h
            for dataset in ("extras","gas")
                HDF5.delete_object(h,dataset)
                h[dataset] = fill(1e-200, 1,1,1,1,1)
            end
        end
        open_snapshot(p) do s
            @test summarize(s,"tracer").mean == 1e-200
            @test_throws DomainError summarize(s,"tracer"; density="rho")
        end
    end
end
