# Run each worker/chunk configuration in a fresh Julia process.
include(joinpath(@__DIR__, "..", "examples", "mhd_workflow.jl"))
using .MHDWorkflowExample: run_mhd_workflow, write_demo

function time_workflow(dir, output, workers, chunk_blocks, phase)
    GC.gc()
    previous_peak = Int(Sys.maxrss())
    trial = @timed run_mhd_workflow(dir, output; workers, chunk_blocks)
    peak = Int(Sys.maxrss())
    r = trial.value
    println((; phase, workers, chunk_blocks, threads=Threads.nthreads(), cells=r["selected_cells"],
             ms=trial.time*1000, million_cells_per_second=r["selected_cells"]/trial.time/1e6,
             allocated_MiB=trial.bytes/2.0^20, peak_before_MiB=previous_peak/2.0^20,
             peak_after_MiB=peak/2.0^20, kinetic=r["kinetic_energy"], magnetic=r["magnetic_energy"]))
end

function main(args)
    length(args) == 3 || error("Usage: julia -t 4 --project=. benchmark/workflow.jl RUN_DIRECTORY WORKERS CHUNK_BLOCKS")
    dir, w, c = args
    workers, chunk = parse(Int,w), parse(Int,c)
    mktempdir() do tmp
        demo = mkdir(joinpath(tmp, "demo"))
        write_demo(demo)
        time_workflow(demo, joinpath(tmp,"synthetic"), workers, chunk, "synthetic warmup")
        time_workflow(dir, joinpath(tmp,"first"), workers, chunk, "first real")
        for i in 1:3
            time_workflow(dir, joinpath(tmp,"trial$i"), workers, chunk, "warm $i")
        end
    end
end

main(ARGS)
