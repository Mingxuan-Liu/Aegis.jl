@testset "Newtonian MHD quantities" begin
    mktempdir() do dir
        names = ("rho","press","vel1","vel2","vel3","Bcc1","Bcc2","Bcc3","tracer")
        p = write_fixture(joinpath(dir,"gas.athdf"); fields=names, blocksize=(3,2,1))
        h5open(p,"r+") do file
            file["gas"][:,:,:,:,1] = fill(4f0,3,2,1,2)
            for (n,value) in enumerate((6,3,4,0,0,0,2,11))
                file["extras"][:,:,:,:,n] = fill(Float32(value),3,2,1,2)
            end
        end
        open_snapshot(p) do s
            # Choose scales so T in K equals P/rho in these manufactured units.
            g = NewtonianMHD(eos=IdealGas(5/3),velocity_unit_cgs=1e5,
                             mean_particle_mass_g=1.380649e-26)
            answers = (speed=5.0, magnetic_magnitude=2.0, kinetic_energy_density=50.0,
                magnetic_energy_density=2.0, magnetic_pressure=2.0, gas_pressure=6.0,
                plasma_beta=3.0, alfven_speed=1.0, thermal_energy_density=9.0,
                sound_speed=sqrt(2.5), temperature=1.5, sonic_mach=sqrt(10.0), alfven_mach=5.0)
            mesh = CartesianMesh(s)
            for (name,answer) in pairs(answers)
                q = quantity(s,g,name)
                slab = read_field(s,q,2; ranges=(2:3,2:2,1:1))
                @test size(slab) == (2,1,1)
                @test eltype(slab) === Float64
                @test all(v -> isapprox(v,answer; rtol=6eps()),slab)
                result = summarize(mesh,q; workers=2)
                @test result.mean ≈ answer rtol=8eps()
                @test result.integral ≈ 63answer rtol=8eps()
            end
            @test quantity(s,g,:temperature).units == "K"
            # An explicitly remapped custom scalar, and unrelated missing inputs.
            custom = quantity(s,g,:gas_pressure; pressure="tracer",density="absent")
            @test only(unique(read_field(s,custom,1))) == 11
            @test_throws ArgumentError quantity(s,:relativistic,:kinetic_energy_density)
            @test_throws ArgumentError quantity(s,NewtonianMHD(),:sound_speed)
            @test_throws ArgumentError quantity(s,NewtonianMHD(eos=IdealGas(5/3)),:temperature)
            @test_throws ArgumentError quantity(s,g,:alfven_speed; density="absent")
            @test_throws ArgumentError quantity(s,g,:unknown)
            iso = NewtonianMHD(eos=Isothermal(2.5),velocity_unit_cgs=1e5,
                               mean_particle_mass_g=1.380649e-26)
            for (name,want) in pairs((gas_pressure=25.,plasma_beta=12.5,sound_speed=2.5,
                                      sonic_mach=2.,temperature=6.25))
                q = quantity(s,iso,name; pressure="missing")
                @test all(x -> isapprox(x,want;rtol=6eps()),read_field(s,q,2))
            end
            @test isempty(quantity(s,iso,:sound_speed).inputs)
            @test_throws ArgumentError read_field(s,quantity(s,iso,:sound_speed),3)
            @test_throws ArgumentError read_field(s,quantity(s,iso,:sound_speed),1; ranges=(1:4,1:1,1:1))
            @test_throws ArgumentError quantity(s,iso,:thermal_energy_density)
        end
        h5open(p,"r+") do file
            file["extras"][:,:,:,:,7] = zeros(Float32,3,2,1,2)
            file["extras"][1,1,1,1,1] = -1f0
            file["gas"][1,1,1,1,1] = 0f0
        end
        open_snapshot(p) do s
            g = NewtonianMHD(eos=IdealGas(1.4))
            @test all(iszero,read_field(s,quantity(s,g,:alfven_speed),2))
            @test all(iszero,read_field(s,quantity(s,g,:magnetic_energy_density),2))
            @test_throws DomainError read_field(s,quantity(s,g,:plasma_beta),2)
            @test_throws DomainError read_field(s,quantity(s,g,:alfven_mach),2)
            @test_throws DomainError read_field(s,quantity(s,g,:kinetic_energy_density),1)
            @test_throws DomainError read_field(s,quantity(s,g,:thermal_energy_density),1)
        end
        # Rebuild in Float64 to exercise primitive scales outside Float32.
        h5open(p,"r+") do file
            gas = fill(1e300,size(file["gas"]))
            extra = Float64.(read(file["extras"]))
            extra[:,:,:,:,1] .= 1e300
            extra[:,:,:,:,7] .= 1e200
            HDF5.delete_object(file,"gas")
            HDF5.delete_object(file,"extras")
            file["gas"] = gas
            file["extras"] = extra
        end
        open_snapshot(p) do s
            # The ratio is finite mathematically. Its intermediate field
            # pressure exceeds Float64 range; a zero beta would be spurious.
            @test_throws DomainError read_field(s,quantity(s,NewtonianMHD(),:plasma_beta),2)
            iso = NewtonianMHD(eos=Isothermal(1e-315))
            @test_throws DomainError read_field(s,quantity(s,iso,:gas_pressure),2)
        end
        h5open(p,"r+") do file
            file["extras"][:,:,:,:,1] = fill(1e-300,3,2,1,2)
        end
        open_snapshot(p) do s
            g = NewtonianMHD(eos=IdealGas(1.4),velocity_unit_cgs=1e150,mean_particle_mass_g=1.380649e-16)
            @test_throws DomainError read_field(s,quantity(s,g,:temperature),2)
        end
        no_rho = write_fixture(joinpath(dir,"single.athdf"); fields=("q",))
        open_snapshot(no_rho) do s
            iso = NewtonianMHD(eos=Isothermal(2.))
            @test all(==(2.),read_field(s,quantity(s,iso,:sound_speed),1))
            @test_throws ArgumentError quantity(s,iso,:kinetic_energy_density)
        end
    end
    @test_throws ArgumentError IdealGas(1)
    @test_throws ArgumentError Isothermal(-1)
    @test_throws ArgumentError NewtonianMHD(eos=:relativistic)
    @test_throws ArgumentError NewtonianMHD(velocity_unit_cgs=NaN)
end

@testset "mixed-precision derived slab allocation" begin
    mktempdir() do dir
        shape = (32,32,32)
        p = write_fixture(joinpath(dir,"block.athdf"); order=[1], blocksize=shape,
                          fields=("vel1","vel2","vel3","rho"))
        # Exercise the runtime input-array types. Density has a small signal
        # below the Float32 precision of the other three components.
        h5open(p,"r+") do h
            density = Float64.(read(h,"gas")) .+ 2.0^-31
            HDF5.delete_object(h,"gas")
            h["gas"] = density
        end
        open_snapshot(p) do s
            kinetic = quantity(s,NewtonianMHD(),:kinetic_energy_density)
            slab = read_field(s,kinetic,1)
            @test first(slab) ≈ (2 + 2.0^-31)*1.5 rtol=4eps(Float64)
            # Warmed end-to-end call: three Float32 velocities, Float64
            # density and Float64 output. Allow 1 MiB additional HDF5.jl
            # bookkeeping without permitting arrays of boxed cell values.
            payload = prod(shape)*(3sizeof(Float32)+2sizeof(Float64))
            alloc = @allocated read_field(s,kinetic,1)
            @test alloc < payload + 2^20
        end
    end
end
