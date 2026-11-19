# Manufacture point samples at prescribed global cell centers. A permutation of
# the entire 3D block lattice exercises lookup across all three coordinate axes.
function derivative_fixture(path, cells, layout, fieldfn; irregular=false)
    @assert all(cells .% layout .== 0)
    ns = cells .÷ layout
    locations = reverse(vec(collect(CartesianIndices(layout))))
    fields = ("q", "u", "v", "w", "const", "B1")
    write_fixture(path; order=collect(1:length(locations)), fields, blocksize=ns)
    faces = ntuple(3) do d
        t = collect(range(0.0, 1.0; length=cells[d]+1))
        irregular ? d .* (t .+ t.^2) ./ 2 : t
    end
    centers = map(f -> (f[1:end-1] .+ f[2:end]) ./ 2, faces)
    h5open(path, "r+") do h
        function resetattr(k,v)
            HDF5.delete_attribute(h,k)
            attributes(h)[k] = v
        end
        resetattr("RootGridSize", Int32[cells...])
        h["LogicalLocations"][:,:] = reduce(hcat, [collect(Tuple(l)) .- 1 for l in locations])
        for d in 1:3
            resetattr("RootGridX$d", [0.,last(faces[d]),irregular ? -1. : 1.])
            for (b,l) in enumerate(locations)
                firstcell = (l[d]-1)*ns[d] + 1
                h["x$(d)f"][:,b] = faces[d][firstcell:firstcell+ns[d]]
                h["x$(d)v"][:,b] = centers[d][firstcell:firstcell+ns[d]-1]
            end
        end
        data = zeros(Float64, ns..., length(locations), length(fields))
        for (b,l) in enumerate(locations), idx in CartesianIndices(ns)
            x = ntuple(d -> centers[d][(l[d]-1)*ns[d]+idx[d]], 3)
            data[idx,b,:] .= (fieldfn(x...)..., 42., 2.)
        end
        HDF5.delete_object(h,"extras")
        h["extras"] = data
    end
    path
end

@testset "quadratic interpolation at centers and across irregular Cartesian blocks" begin
    mktempdir() do dir
        function affine(x,y,z)
            (x^2 + 2y^2 + 4x*z-y*z, x+2y+3z, 4x+5y+6z, 7x+8y+9z)
        end
        p = derivative_fixture(joinpath(dir,"irregular.athdf"),(6,4,3),(3,2,3),affine;irregular=true)
        open_snapshot(p) do s
            mesh = CartesianMesh(s)
            grid = DifferentialGrid(mesh; boundary=(:onesided,:onesided,:onesided))
            tensor = vector_gradient(grid, ("u","v","w"))
            vort = curl(grid, ("u","v","w"))
            g = gradient(grid, interpret_field(s,"q";units="K"))
            div = divergence(grid, ("u","v","w"))
            @test g[1].units == "(K)/code_length"
            # Absolute errors scale with spacing and roundoff in the stored
            # linear combinations; all analytic components have magnitude O(10).
            for b in 1:s.nblocks
                x,y,z = block_coordinates(s,b).centers
                expected = ([2a+4c for a=x,b=y,c=z], [4b-c for a=x,b=y,c=z],
                            [4a-b for a=x,b=y,c=z])
                for j in 1:3
                    @test read_field(s,g[j],b) ≈ expected[j] atol=2e-12
                    @test all(v -> isapprox(v, (2.,-4.,2.)[j]; atol=2e-12),
                              read_field(s,vort[j],b))
                    for i in 1:3
                        @test all(v -> isapprox(v, 3(i-1)+j; atol=2e-12),
                                  read_field(s,tensor[i][j],b))
                    end
                end
                @test all(v -> isapprox(v,15.;atol=2e-12), read_field(s,div,b))
                @test all(iszero,read_field(s,gradient(grid,"const")[1],b))
            end
            @test summarize(mesh,div).mean ≈ 15
            @test_throws ArgumentError DifferentialGrid(mesh)
            @test_throws ArgumentError DifferentialGrid(mesh;boundary=(:invariant,:onesided,:onesided))
            @test_throws ArgumentError gradient(grid,"B1")
            @test_throws ArgumentError divergence(grid, (interpret_field(s,"u";units="velocity"),
                "v",interpret_field(s,"w";units="pressure")))
            @test_throws ArgumentError read_field(s,g[1],1;ranges=(0:1,1:1,1:1))
        end
    end
end

@testset "periodic curl convergence by location and block partition" begin
    mktempdir() do dir
        errors = zeros(3,3)
        A = 2.25
        k = 2π
        wave(x,y,z) = (0.,0.,A*sin(k*x),0.)
        for (res,n) in enumerate((24,48,96))
            partitions = Vector{Float64}[]
            for nb in (1,4)
                p = derivative_fixture(joinpath(dir,"wave.athdf"), (n,1,1),(nb,1,1),wave)
                recovered = zeros(n)
                open_snapshot(p) do s
                    mesh = CartesianMesh(s)
                    grid = DifferentialGrid(mesh;boundary=(:periodic,:invariant,:invariant))
                    ω = curl(grid, ("u","v","w"))
                    for b in 1:s.nblocks
                        idx = s.locations[1,b]*s.blocksize[1] .+ (1:s.blocksize[1])
                        recovered[idx] = vec(read_field(s,ω[3],b))
                    end
                    # Select just the two cells straddling a seam. Their outer
                    # neighbors fall outside the requested interval.
                    window=Region(x1=(0.5-1/n,0.5+1/n))
                    selected = Float64[]
                    for b in selectedblocks(mesh;region=window)
                        append!(selected,vec(read_field(s,ω[3],b.block;ranges=b.ranges)))
                    end
                    @test sort(selected) == sort(recovered[n÷2:n÷2+1])
                    @test all(iszero,read_field(s,gradient(grid,"v")[3],1))
                    @test_throws ArgumentError DifferentialGrid(mesh; boundary=(:periodic,:periodic,:periodic))
                end
                push!(partitions,recovered)
            end
            @test only(unique(partitions)) == partitions[1]
            dx = 1/n
            pointwise = abs.(partitions[1] .- A*k.*cos.(k .* ((1:n) .- 0.5).*dx))
            wall = [1,n]
            seams = [i for i in 2:n-1 if mod(i,n÷4) in (0,1)]
            interior = setdiff(1:n,wall,seams)
            for (j,category) in enumerate((interior,seams,wall))
                errors[res,j] = maximum(pointwise[category])
            end
        end
        rates = log2.(errors[1:end-1,:] ./ errors[2:end,:])
        @test all(r -> 1.75 < r < 2.15, rates)
    end
end

@testset "short active axis and invalid halo data" begin
    mktempdir() do dir
        f(x,y,z)=(x+y+z,0.,0.,0.)
        p=derivative_fixture(joinpath(dir,"short.athdf"),(2,1,1),(2,1,1),f)
        open_snapshot(p) do s
            @test_throws ArgumentError DifferentialGrid(CartesianMesh(s);boundary=(:periodic,:invariant,:invariant))
        end
        p=derivative_fixture(joinpath(dir,"nan.athdf"),(6,1,1),(3,1,1),f)
        h5open(p,"r+") do h
            # Last sample of the middle block is the left neighbor in file block 1.
            h["extras"][2,1,1,2,1] = NaN
        end
        open_snapshot(p) do s
            grid = DifferentialGrid(CartesianMesh(s);boundary=(:periodic,:invariant,:invariant))
            @test_throws DomainError read_field(s,gradient(grid,"q")[1],1;ranges=(1:1,1:1,1:1))
        end
    end
end
