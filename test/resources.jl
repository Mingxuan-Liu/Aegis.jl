# A small hyperslab must stay cheaper than materializing its containing field.
# Warm the identical call first so compiler allocations are outside the check.
_read_probe(s, ranges) = read_field(s, "tracer", 1; ranges)

@testset "selected read allocation" begin
    mktempdir() do dir
        shape = (128, 96, 32)
        path = write_fixture(joinpath(dir, "large-block.athdf"); order=[1],
                             blocksize=shape, fields=("tracer",))
        open_snapshot(path) do s
            ranges = (100:101, 2:3, 1:2)
            @test size(_read_probe(s, ranges)) == (2, 2, 2)
            allocated = @allocated _read_probe(s, ranges)
            # Leave ample space for HDF5.jl bookkeeping across supported Julia
            # versions. The whole field is 1.5 MiB; its eight-cell slab is 32 B.
            @test allocated < prod(shape) * sizeof(Float32) ÷ 4
        end
    end
end
