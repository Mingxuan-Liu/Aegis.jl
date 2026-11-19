module WorkflowTests

using Test
using Aegis
using HDF5
using TOML
include(joinpath(@__DIR__, "..", "examples", "mhd_workflow.jl"))
using .MHDWorkflowExample: run_mhd_workflow, write_demo

@testset "Newtonian workflow and persistent products" begin
    mktempdir() do tmp
        run = mkdir(joinpath(tmp, "input"))
        write_demo(run)
        for (i,(workers,chunk)) in enumerate(((1,1),(2,1),(2,3)))
            out = joinpath(tmp, "out$i")
            m = run_mhd_workflow(run,out; workers,chunk_blocks=chunk,slice_field="tracer")
            @test m["selected_cells"] == 24
            @test m["selected_volume"] == 189
            @test m["kinetic_energy"] == 12600
            @test m["magnetic_energy"] == 378
            @test m["included_mass"] == 1008
            @test m["excluded_mass"] == 0
            @test TOML.parsefile(joinpath(out,"summary.toml"))["probe_value"] == m["probe_value"]
            h5open(joinpath(out,"mass_by_density.h5"),"r") do h
                @test sum(read(h,"totals")) ≈ 1008
            end
            h5open(joinpath(out,"slice.h5"),"r") do h
                @test read(h,"values") == repeat(collect(1:6),1,2)
            end
        end
        rm(joinpath(run,"flow.hst"))
        rm(joinpath(run,"athinput.flow"))
        # The last three x1 cells are wider than the first three. Restrict to
        # the first block in physical space, which is last in file order.
        m = run_mhd_workflow(run,joinpath(tmp,"meshonly"); region=Region(x1=(0,6)))
        @test isempty(m["histories"]) && isempty(m["inputs"])
        @test m["selected_volume"] == 54
        @test m["included_mass"] == 180
        @test m["kinetic_energy"] == 2250
        @test m["magnetic_energy"] == 108
        @test m["slice_field"] == "magnetic_energy_density"
        narrow = run_mhd_workflow(run,joinpath(tmp,"narrow"); slice_axis=1,
                                  slice_field=s->quantity(s,NewtonianMHD(),:speed), region=Region(x1=(0.49,2)))
        @test narrow["slice_at"] == 0.5
        @test narrow["probe_value"] == 5.0
        @test_throws ArgumentError run_mhd_workflow(run,joinpath(tmp,"meshonly"))
        failed = joinpath(tmp,"failed")
        @test_throws ArgumentError run_mhd_workflow(run,failed; slice_field="missing")
        @test !ispath(failed)
    end
end
end
