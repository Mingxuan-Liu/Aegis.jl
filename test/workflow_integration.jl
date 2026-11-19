# Private workflow validation. All three supplied runs are required.
using Test
using HDF5
using Aegis
using TOML
include(joinpath(@__DIR__, "..", "examples", "mhd_workflow.jl"))
using .MHDWorkflowExample: run_mhd_workflow, last_snapshot

function check_history_energy(h, out, key, suffix)
    times = history_column(h, "time")
    row = only(findall(==(out["time"]), times))
    # Each direction is an independent six-significant-digit printed integral.
    components = [history_column(h, "$d-$suffix")[row] for d in 1:3]
    atol = sum(x -> iszero(x) ? 0.0 : 0.5*10.0^(floor(log10(abs(x)))-5), components)
    # Nonrelativistic rho, velocity and cell-centered B are saved at Float32.
    # Up to three independently rounded factors enter a kinetic product.
    atol += 4eps(Float32)*sum(abs,components)
    @test abs(out[key] - sum(components)) <= atol
    sum(components)
end

function workflow_checks(parent)
    for (name,is_pic) in (("3D_Athena/128_eta5e-3",false), ("3D_Athena/256_eta5e-3",false),
                          ("2D_MHD_PIC/run424",true))
        dir = joinpath(parent,name)
        isdir(dir) || error("Missing required workflow integration run: $dir")
        _, ref = last_snapshot(dir, nothing)
        # Interpret selected named custom output as an uninterpreted array here.
        # Both MHD outputs have a named heating field. The PIC mesh stores fluid
        # primitives only; a magnetic-energy slice exercises its compatible model.
        slice_field = open_snapshot(ref.path) do s
            is_pic && return nothing
            only([f.name for f in s.fields if f.name == "ohmic_heating"])
        end
        mktempdir() do tmp
            path = joinpath(tmp,"workflow")
            r = run_mhd_workflow(dir, path; slice_field, workers=min(4,Threads.nthreads()), chunk_blocks=2)
            h = read_history(only(r["histories"])["path"])
            ke = check_history_energy(h,r,"kinetic_energy","KE")
            me = check_history_energy(h,r,"magnetic_energy","ME")
            mass_row = only(findall(==(r["time"]), history_column(h,"time")))
            mass = history_column(h,"mass")[mass_row]
            @test isapprox(r["included_mass"] + r["excluded_mass"],mass; rtol=5e-6)
            stored = TOML.parsefile(joinpath(path,"summary.toml"))
            @test stored["source"] == ref.path
            @test stored["selected_cells"] == prod(stored["root_cells"])
            @test stored["magnetic_energy"] == r["magnetic_energy"]
            for filename in ("mass_by_density.h5","slice.h5")
                @test isfile(joinpath(path,filename))
            end
            println((; name, t=r["time"], kinetic=r["kinetic_energy"], solver_kinetic=ke,
                     magnetic=r["magnetic_energy"], solver_magnetic=me, slice=r["slice_field"],
                     included_mass=r["included_mass"], excluded_mass=r["excluded_mass"]))
        end
    end
end

length(ARGS) <= 1 || error("usage: test/workflow_integration.jl [private-data-root]")
@testset "Full workflow on private fluid outputs" workflow_checks(isempty(ARGS) ? joinpath(@__DIR__,"test_data") : only(ARGS))
