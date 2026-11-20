@testset "Magnetic fluctuation energy about a declared guide" begin
    mktempdir() do dir
        path = joinpath(dir,"wave.athdf")
        order = [2,1]
        write_fixture(path; order, fields=("Bx","By","Bz","uperp"), blocksize=(4,2,1))
        # A physical circular wave on shuffled, uneven blocks. Rotate all
        # three axes together by a proper orthogonal matrix with rational entries.
        rotate((x,y,z)) = ((x-2y-2z)/3, (-2x+y-2z)/3, (2x+2y-z)/3)
        for rotation in (identity, rotate)
            h5open(path, "r+") do h
                x = read(h["x1v"])
                values = zeros(4,2,1,2,4)
                for b in 1:2, i in 1:4
                    phase = 2π*x[i,b]/36
                    for (slot,v) in pairs(rotation((1.25, 3cos(phase), 3sin(phase))))
                        values[i,:,1,b,slot] .= v
                    end
                end
                values[:,:,:,:,4] .= 4.5
                HDF5.delete_object(h,"extras")
                h["extras"] = values
            end
            open_snapshot(path) do s
                gas = NewtonianMHD()
                reference = rotation((1,0,0))
                direction = rotation((2,0,0))
                magnetic = ("Bx", "By", "Bz")
                fields = map((:parallel, :transverse)) do part
                    magnetic_fluctuation(s,gas,part; reference,direction,magnetic)
                end
                mesh = CartesianMesh(s)
                # Parallel perturbation 1/4, transverse amplitude 3. The
                # 48-eps relative margin on the small parallel energy accounts
                # for constructing the field in a rotated Float64 basis.
                for (f, expected) in zip(fields, (1/32, 4.5))
                    for b in 1:2
                        slab = read_field(s,f,b; ranges=(2:4,1:1,1:1))
                        @test size(slab) == (3,1,1)
                        @test all(e -> isapprox(e,expected; rtol=48eps()),slab)
                    end
                    @test summarize(mesh,f).integral ≈ expected*108 rtol=48eps()
                end
                ref = interpret_field(s,"uperp"; description="in-situ wave energy")
                stored = magnetic_fluctuation(s,gas,:transverse; reference,direction,stored=ref)
                cut = slice_plane(mesh,stored; axis=2,at=1.5)
                @test all(==(4.5), cut.values)
                @test cut.metadata.field.source === :stored
                @test only(cut.metadata.field.inputs).description == "in-situ wave energy"
                @test only(cut.metadata.field.inputs).dataset == "extras"
                calc = slice_plane(mesh,fields[2]; axis=2,at=1.5)
                @test calc.metadata.field.source === :reconstructed
                @test calc.metadata.field.inputs[2].name == "By"
                @test calc.metadata.field.reference == reference
                @test calc.metadata.field.direction == direction
                @test calc.metadata.field.frame == :simulation
                @test calc.metadata.field.permeability == 1
                @test_throws ArgumentError magnetic_fluctuation(s,gas,:transverse; reference,direction,magnetic,stored=ref)
                @test_throws ArgumentError magnetic_fluctuation(s,gas,:transverse; reference=(0,0,0),stored=ref)
                @test_throws ArgumentError magnetic_fluctuation(s,gas,:parallel; reference,direction=(1,NaN,0),magnetic)
                @test_throws ArgumentError magnetic_fluctuation(s,gas,:parallel; reference,magnetic=("Bx","By","By"))
                @test_throws ArgumentError magnetic_fluctuation(s,gas,:parallel; reference,magnetic=("Bx","By",stored))
                @test_throws ArgumentError magnetic_fluctuation(s,gas,:helicity; reference,magnetic)
            end
        end

        # Total-minus-parallel magnetic energy loses the entire transverse
        # perturbation here. The x-y difference gives an exact independent oracle.
        h5open(path, "r+") do h
            field = h["extras"]
            field[1,1,1,1,1] = 2.0^45
            field[1,1,1,1,2] = 2.0^45+1
            field[1,1,1,1,3] = 0.0
            field[1,1,1,1,4] = -1.0
            field[2,1,1,1,1] = Inf
            field[3,1,1,1,1] = 2.0^54
            field[3,1,1,1,2] = 2.0^54
            field[3,1,1,1,3] = 0
        end
        open_snapshot(path) do s
            gas = NewtonianMHD()
            opts = (; reference=(0,0,0), magnetic=("Bx","By","Bz"))
            for d in (1.0, floatmax(Float64), nextfloat(0.0))
                q = magnetic_fluctuation(s,gas,:transverse; opts..., direction=(d,d,0))
                @test only(read_field(s,q,1; ranges=(1:1,1:1,1:1))) ≈ 1/4 rtol=3eps()
                @test_throws DomainError read_field(s,q,1; ranges=(2:2,1:1,1:1))
            end
            # The exact difference between the large Cartesian components
            # is -1 after reference subtraction, even though 2^54-1 rounds.
            tip = magnetic_fluctuation(s,gas,:parallel; reference=(1,0,0),
                                       direction=(1,-1,0), magnetic=opts.magnetic)
            @test only(read_field(s,tip,1; ranges=(3:3,1:1,1:1))) ≈ 1/4 rtol=3eps()
            q = magnetic_fluctuation(s,gas,:parallel; reference=(1,0,0), stored="uperp")
            @test_throws DomainError read_field(s,q,1; ranges=(1:1,1:1,1:1))
        end
    end
end

@testset "Bounded guide-field slab allocation" begin
    mktempdir() do dir
        shape = (32,32,32)
        p = write_fixture(joinpath(dir,"guide.athdf"); order=[1], blocksize=shape,
                          fields=("Bcc1","Bcc2","Bcc3"))
        open_snapshot(p) do s
            field = magnetic_fluctuation(s,NewtonianMHD(),:transverse; reference=(0.25,-0.125,0.0625),
                                         direction=(3,4,0))
            read_field(s,field,1)
            # One Float64 output and three Float32 input arrays; allow a
            # further MiB of HDF5/descriptor overhead, excluding per-cell boxing.
            payload = prod(shape)*(3sizeof(Float32) + sizeof(Float64))
            @test (@allocated read_field(s,field,1)) < payload + 2^20
        end
    end
end
