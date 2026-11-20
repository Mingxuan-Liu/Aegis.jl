using Aegis, Test

@testset "pinned Athena++ SRMHD writer, uniform oblique field" begin
    open_snapshot(joinpath(@__DIR__, "oblique.prim.00000.athdf")) do prim
        open_snapshot(joinpath(@__DIR__, "oblique.cons.00000.athdf")) do cons
            @test prim.time == cons.time == 0.0
            @test prim.cycle == cons.cycle == 0
            @test prim.locations == cons.locations
            @test prim.nblocks == cons.nblocks == 2
            @test prim.blocksize == cons.blocksize == (4, 1, 1)
            gas = SpecialRelativisticMHD()
            # The generator accepts v, while the writer labels u = Wv as vel1:3.
            primitives = ("rho" => 4.0, "press" => 1.0, "vel1" => 0.45,
                          "vel2" => 0.6, "vel3" => 0.0, "Bcc1" => 3.0,
                          "Bcc2" => 4.0, "Bcc3" => 5.0)
            diagnostics = (:lorentz_factor => 1.25, :lab_velocity_1 => 0.36,
                           :lab_velocity_2 => 0.48, :lab_velocity_3 => 0.0,
                           :comoving_magnetic_pressure => 20.5)
            # The initial writer includes the solver's conserved-to-primitive recovery.
            # Allow a few Float64 ulps, including an absolute scale for the zero component.
            agrees(x, want) = isapprox(x, want; rtol=16eps(), atol=2eps())
            for block in 1:prim.nblocks
                @test block_coordinates(prim, block) == block_coordinates(cons, block)
                for (field, expected) in primitives
                    @test all(x -> agrees(x, expected), read_field(prim, field, block))
                end
                for (name, expected) in diagnostics
                    q = quantity(prim, gas, name)
                    @test all(x -> agrees(x, expected), read_field(prim, q, block))
                end
                D = read_field(cons, "dens", block)
                @test all(==(5.0), D)
                W = read_field(prim, quantity(prim, gas, :lorentz_factor), block)
                rho = read_field(prim, "rho", block)
                @test all(agrees.(rho .* W, D))
            end
        end
    end
end
