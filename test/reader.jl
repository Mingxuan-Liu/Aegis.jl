@testset "native blocks and original field identities" begin
    mktempdir() do dir
        p = write_fixture(joinpath(dir, "grid.athdf"); fields=("tracer","marker","marker","B2"))
        h5open(p, "r+") do file
            attributes(file)["run_tag"] = "manufactured"
            # Two slots may have the same advertised name and different data.
            file["extras"][:,:,:,1,2] = fill(-10f0, 3,2,2)
        end
        escaped, data = open_snapshot(p) do s
            @test s.time == 1.25
            @test s.cycle == 8
            @test s.attributes["run_tag"] == "manufactured"
            @test s.blocksize == (3,2,2)
            @test s.locations[:,1] == [1,0,0]
            @test s.levels == [0,0]
            @test [f.name for f in s.fields] == ["tracer","marker","marker","B2"]
            @test_throws ArgumentError resolve_field(s, "rho")
            @test_throws ArgumentError resolve_field(s, "marker")
            @test_throws ArgumentError resolve_field(s, ("extras", 0))
            @test resolve_field(s, ("extras", 2)).name == "marker"
            @test all(==(-10f0), read_field(s, ("extras", 2), 1))
            @test resolve_field(s, "B2").centering == :unknown
            coords = block_coordinates(s, 1)
            @test coords.faces == ([6.,10.,15.,21.], [0.,1.,3.], [0.,1.,3.])
            @test coords.centers[1] == [8.,12.5,18.]
            slab = read_field(s, s.fields[3], 1; ranges=(2:3,1:2,2:2))
            @test eltype(slab) == Float32
            @test size(slab) == (2,2,1)
            narrow = (Int32(2):Int32(3), UInt(1):UInt(2), big(2):big(2))
            @test read_field(s, s.fields[3], Int32(1); ranges=narrow) == slab
            @test_throws ArgumentError read_field(s, "tracer", 1;
                ranges=(big(1):big(10)^50, 1:2, 1:2))
            @test slab[:,:,1] == Float32[215 225; 216 226]
            @test read_field(s,"tracer",2)[:,1,1] == Float32[1,2,3]
            for ranges in ((0:1,1:1,1:1), (1:1,1:3,1:1), (2:1,1:1,1:1), (1:1:3,1:2:2,1:1))
                @test_throws ArgumentError read_field(s, "tracer", 1; ranges)
            end
            @test_throws ArgumentError read_field(s, "tracer", 0)
            @test_throws ArgumentError block_coordinates(s, 3)
            s, slab
        end
        @test data[1] == 215
        @test_throws ArgumentError read_field(escaped, "tracer", 1)
        open_snapshot(write_fixture(joinpath(dir, "other.athdf"))) do s
            @test_throws ArgumentError read_field(s, escaped.fields[1], 1)
        end

        p2 = write_fixture(joinpath(dir, "double.athdf"); fields=("marker",), blocksize=(4,3,1))
        h5open(p2, "r+") do f
            doubles = Float64.(read(f, "extras")) .+ 2.0^-44
            HDF5.delete_object(f, "extras")
            f["extras"] = doubles
        end
        open_snapshot(p2) do s
            a = read_field(s, "marker", 2)
            @test eltype(a) == Float64
            @test size(a) == (4,3,1)
            @test a[4,3,1] == 134.0 + 2.0^-44
        end
    end
end

@testset "malformed native metadata" begin
    mktempdir() do dir
        function altered(edit!)
            p = write_fixture(joinpath(dir, "broken.athdf"))
            h5open(edit!, p, "r+")
            open_snapshot(identity, p)
        end
        @test_throws ArgumentError altered(f -> HDF5.delete_attribute(f, "Time"))
        @test_throws ArgumentError altered(f -> HDF5.delete_object(f, "x2v"))
        @test_throws ArgumentError altered() do f
            HDF5.delete_attribute(f, "NumVariables")
            attributes(f)["NumVariables"] = Int32[2,3]
        end
        @test_throws ArgumentError altered() do f
            HDF5.delete_object(f, "extras")
            f["extras"] = zeros(Float32, 2,3,2,2,2) # swapped x1 and x2
        end
        @test_throws ArgumentError altered() do f
            HDF5.delete_object(f, "extras")
            f["extras"] = zeros(UInt16, 3,2,2,2,2) # packed/quantized data
        end
        # The time is sufficient to discover a stream with no recorded cycle.
        @test altered(f -> HDF5.delete_attribute(f, "NumCycles")).cycle === nothing
    end
end

@testset "discover separate streams by validated time" begin
    mktempdir() do dir
        for (file, t) in (("run.out2.00000.athdf",3.), ("run.out2.00002.athdf",1.),
                          ("run.out2.00001.athdf",1.), ("run.out1.00001.athdf",9.),
                          ("extra.athdf",9.))
            write_fixture(joinpath(dir, file); time=t)
        end
        streams = discover_snapshots(dir)
        @test Set(keys(streams)) == Set(["run.out1","run.out2","extra.athdf"])
        @test [r.time for r in streams["run.out2"]] == [1.,1.,3.]
        @test [basename(r.path) for r in streams["run.out2"]] ==
            ["run.out2.00001.athdf","run.out2.00002.athdf","run.out2.00000.athdf"]
        broken = joinpath(dir, "broken.athdf")
        write_fixture(broken; time=NaN)
        failure = try
            discover_snapshots(dir)
        catch e
            e
        end
        @test failure isa ArgumentError
        @test occursin(broken, sprint(showerror, failure))
        @test occursin("Time", sprint(showerror, failure))
        rm(broken)
        write_fixture(joinpath(dir, "extra.athdf.00000.athdf"))
        @test_throws ArgumentError discover_snapshots(dir)
    end
end
