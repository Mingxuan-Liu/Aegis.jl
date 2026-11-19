# Run from the repository root. See docs/performance.md for the measurement protocol.
using Aegis
using HDF5
include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

function measure(label, f)
    GC.gc()
    before = Int(Sys.maxrss())
    timing = @timed f()
    after = Int(Sys.maxrss())
    println((; phase=label, seconds=timing.time, allocated_bytes=timing.bytes,
             peak_before_bytes=before, peak_after_bytes=after))
    return timing.value
end

# A separately timed I/O comparison. Consume every slab, but perform no reduction.
function readpass(s, field)
    last_value = 0.0
    for b in 1:s.nblocks
        last_value = last(read_field(s, field, b))
    end
    return last_value
end

function workload(s, field, workers, chunk, phase)
    if phase == "summary"
        return summarize(s, field; workers, chunk_blocks=chunk).integral
    elseif phase == "read"
        return readpass(s, field)
    elseif phase == "partial"
        r = ntuple(d -> 1:min(2, s.blocksize[d]), 3)
        return sum(read_field(s, field, 1; ranges=r))
    elseif phase == "mesh"
        CartesianMesh(s)
        return s.nblocks
    end
    error("phase must be summary, read, partial or mesh")
end

# Feed one preloaded block through the production kernel repeatedly. A separate
# scalar accumulator per pass reproduces the serial per-block calculation.
function hot_kernel(block, values, repetitions)
    checksum = 0.0
    for _ in 1:repetitions
        total = Aegis._accumulate!(Aegis.CellTotals(), block, values, nothing)
        checksum += total.integral + total.correction
    end
    checksum
end

function benchmark_kernel(path, fieldname)
    open_snapshot(path) do s
        block = first(selectedblocks(CartesianMesh(s)))
        values = read_field(s, fieldname, block.block; ranges=block.ranges)
        repetitions = 32
        println((; cells_per_trial=length(values)*repetitions))
        hot_kernel(block, values, 1)
        for trial in 1:3
            checksum = measure("steady in-memory kernel trial $trial",
                               () -> hot_kernel(block, values, repetitions))
            println((; checksum))
        end
    end
end

function benchmark(path, fieldname, workers, chunk, phase)
    println((; julia=VERSION, threads=Threads.nthreads(), hdf5=HDF5.API.h5_get_libversion(),
             machine=Sys.MACHINE, snapshot=basename(path), file_bytes=filesize(path),
             field=fieldname, workers, chunk, phase))
    phase == "kernel" && return benchmark_kernel(path, fieldname)
    mktempdir() do tmp
        demo = example_run(tmp)
        measure("warmup (includes compilation)", () -> open_snapshot(demo.snapshot) do s
            workload(s, resolve_field(s, "rho"), workers, chunk, phase)
        end)
    end
    for repeat in 0:3
        label = repeat == 0 ? "first real file (includes format-specific compilation)" : "steady trial $repeat"
        value = measure("$label (open + $phase + close)", () -> open_snapshot(path) do s
            workload(s, resolve_field(s, fieldname), workers, chunk, phase)
        end)
        println((; value))
    end
end

length(ARGS) == 5 || error("Usage: julia --threads=4 --project=. benchmark/native.jl SNAPSHOT FIELD WORKERS CHUNK_BLOCKS PHASE")
benchmark(ARGS[1], ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]), ARGS[5])
