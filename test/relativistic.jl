# Manufactured inertial Cartesian primitives, c=μ=1. These independent states
# use known LAB THREE-velocities, each paired with the corresponding writer u.
# No conserved field in this synthetic file is used as solver evidence.
function sr_test_state(u, v, B, W, pm; rho=4., press=0.1)
    (; u, v, B, W, pm, rho, press)
end

function write_sr_states(path, states, T)
    fields = ("rho", "press", "vel1", "vel2", "vel3", "Bcc1", "Bcc2", "Bcc3", "pB_fluid")
    nx = length(states) ÷ 2
    order = [2, 1]
    @assert 2nx == length(states)
    write_fixture(path; fields, blocksize=(nx,1,1), order)
    h5open(path,"r+") do h
        gas = Array{T}(undef,nx,1,1,2,1)
        extras = Array{T}(undef,nx,1,1,2,8)
        for (block, loc) in enumerate(order), i in 1:nx
            a = states[(loc-1)*nx+i]
            gas[i,1,1,block,1] = a.rho
            extras[i,1,1,block,:] .= (a.press, a.u..., a.B..., a.pm)
        end
        for (name, data) in (("gas",gas), ("extras",extras))
            HDF5.delete_object(h,name)
            h[name] = data
        end
    end
    path
end

@testset "SR ideal-fluid boosts through the primitive reader" begin
    W4 = 1000/sqrt(1999)
    rest = (0.,0.,0.)
    vaxis = (0.6,0.,0.)
    uaxis = (0.75,0.,0.)
    states = [
        sr_test_state(rest,rest,(3.,4.,0.),1.,12.5),
        sr_test_state(uaxis,vaxis,(5.,0.,0.),1.25,12.5),
        sr_test_state(uaxis,vaxis,(0.,5.,0.),1.25,8.),
        sr_test_state((0.45,0.6,0.),(.36,.48,0.),(3.,4.,5.),1.25,20.5),
        # Mignone & Bodo (2006) Table 1, case 4L. Convert the published v to u.
        sr_test_state((0.999W4,0.,0.),(.999,0.,0.),(10.,7.,7.),W4,50.097951; rho=1.),
        sr_test_state((0.,0.,-4/3),(0.,0.,-0.8),(0.,3.,4.),5/3,9.62),
        sr_test_state((0.,0.,-4/3),(0.,0.,-0.8),rest,5/3,0.),
        # Coordinate-permuted and sign-reversed oblique case
        sr_test_state((-.6,0.,-.45),(-.48,0.,-.36),(4.,5.,3.),1.25,20.5)
    ]
    model = SpecialRelativisticMHD()
    vnames = (:lab_velocity_1, :lab_velocity_2, :lab_velocity_3)
    mktempdir() do dir
        for T in (Float32,Float64)
            path = write_sr_states(joinpath(dir,"states-$T.athdf"),states,T)
            open_snapshot(path) do s
                for block in (2,1) # increasing physical x1 despite reversed file order
                    W = vec(read_field(s,quantity(s,model,:lorentz_factor),block))
                    v = map(name -> vec(read_field(s,quantity(s,model,name),block)),vnames)
                    pm = vec(read_field(s,quantity(s,model,:comoving_magnetic_pressure),block))
                    rho = vec(read_field(s,"rho",block))
                    @test eltype(pm) === Float64
                    for i in eachindex(W)
                        a = states[(2-block)*length(W)+i]
                        # A few arithmetic operations and primitive quantization.
                        tol = 12eps(T)
                        @test W[i] ≈ a.W rtol=tol
                        @test all(isapprox(v[d][i],a.v[d];rtol=tol) for d in 1:3)
                        @test pm[i] ≈ a.pm rtol=tol
                        @test rho[i]*W[i] ≈ a.rho*a.W rtol=tol

                        # Boost the fluid-frame magnetic vector (0,B') into
                        # the lab, using the known three-velocity. The ideal
                        # field transforms as B'∥ = B∥ and B'⊥ = B⊥/W.
                        speed = sqrt(sum(abs2,a.v))
                        n = iszero(speed) ? (1.,0.,0.) : a.v ./ speed
                        Bpar = sum(a.B .* n)
                        Bperp = a.B .- Bpar .* n
                        btime = a.W*speed*Bpar
                        bspace = (a.W*Bpar) .* n .+ Bperp ./ a.W
                        usp = ntuple(d -> W[i]*v[d][i],3)
                        # Absolute error grows with cancelling terms in the
                        # Minkowski metric (-,+,+,+), especially for case 4L.
                        @test abs(-W[i]^2 + sum(abs2,usp) + 1) <= tol*W[i]^2
                        @test abs(-W[i]*btime + sum(usp .* bspace)) <= tol*max(1,W[i]*abs(btime))
                        @test abs(2pm[i] - (sum(abs2,bspace)-btime^2)) <= tol*max(1,sum(abs2,bspace))
                    end
                end
            end
        end
        # Binary scaling keeps u·B exactly zero in the written primitives.
        # The transverse pressure remains representable even though separate
        # componentwise rounding of u/W can swamp its small physical value.
        a = 2.0^70
        w = sqrt(1 + 59a^2)
        transverse = sr_test_state((a,3a,7a),(a,3a,7a)./w,(3.,-1.,0.),w,5/w^2)
        extreme = write_sr_states(joinpath(dir,"transverse.athdf"),[transverse,transverse],Float64)
        open_snapshot(extreme) do s
            q = quantity(s,model,:comoving_magnetic_pressure)
            @test only(read_field(s,q,1)) ≈ transverse.pm rtol=16eps()
        end

        # Partial reads on reversed, nonuniform native blocks; selected rows
        # distinguish spatial four-velocity from the returned three-velocity
        # (0.75 vs 0.6) and comoving from lab pressure (8 vs 12.5).
        open_snapshot(joinpath(dir,"states-Float64.athdf")) do s
            pick(q) = only(read_field(s,q,2;ranges=(3:3,1:1,1:1)))
            @test pick("vel1") == .75
            @test pick(quantity(s,model,:lab_velocity_1)) ≈ .6
            @test pick("Bcc2")^2/2 == 12.5
            @test pick(quantity(s,model,:comoving_magnetic_pressure)) ≈ 8
            @test pick(quantity(s,model,:comoving_magnetic_pressure; magnetic=("Bcc2","Bcc1","Bcc3"))) ≈ 12.5
            pfield = interpret_field(s,"pB_fluid"; alias="pB",units="p_unit",description="ideal comoving")
            q = quantity(s,model,:comoving_magnetic_pressure; stored=pfield)
            @test pick(q) == 8
            plane = slice_plane(CartesianMesh(s),q;axis=3,at=0.)
            saved = plane.metadata.field
            @test saved.inputs[1].alias == "pB"
            @test saved.inputs[1].units == "p_unit"
            @test saved.inputs[1].description == "ideal comoving"
            dest = export_product(joinpath(dir,"measured.h5"),plane)
            h5open(dest,"r") do h
                @test read(h["metadata/field/frame"]) == "fluid_rest"
                @test read(h["metadata/field/inputs/1/name"]) == "pB_fluid"
            end
            @test saved.kind === :special_relativistic && saved.source === :stored
            reconstructed = quantity(s,model,:comoving_magnetic_pressure)
            meta = slice_plane(CartesianMesh(s),reconstructed;axis=3,at=0.).metadata.field
            @test meta.source === :reconstructed && length(meta.inputs) == 6
        end
    end
end

@testset "cold slow-flow limit and invalid SR inputs" begin
    # For a transverse weak field, (p_B-p'_B)/p_B = v² exactly.
    # The remainder W-1-v²/2 is (3/8)v⁴+O(v⁶). Use velocities large
    # enough for the quartic remainder to survive Float64 cancellation.
    speeds = (0.,1/25,1/50,1/100)
    B = (0.,1/64,0.)
    pB = sum(abs2,B)/2
    states = map(speeds) do v
        W = inv(sqrt((1-v)*(1+v)))
        sr_test_state((W*v,0.,0.),(v,0.,0.),B,W,pB*(1-v^2);press=1e-8)
    end
    mktempdir() do dir
        model = SpecialRelativisticMHD()
        path = write_sr_states(joinpath(dir,"cold.athdf"),states,Float64)
        open_snapshot(path) do s
            function global_values(name)
                q = quantity(s,model,name)
                vcat(vec(read_field(s,q,2)), vec(read_field(s,q,1)))
            end
            W = global_values(:lorentz_factor)
            pm = global_values(:comoving_magnetic_pressure)
            @test W[1] == 1 && pm[1] ≈ pB
            for i in 2:4
                v = speeds[i]
                @test isapprox((pB-pm[i])/pB,v^2; atol=8eps())
                # The next coefficient in the binomial expansion is 5/16.
                @test abs((W[i]-1-v^2/2)/v^4 - 3/8) < 0.32v^2 + 4eps()/v^4
            end
            @test_throws ArgumentError quantity(s,model,:speed)
            @test_throws ArgumentError quantity(s,model,:lorentz_factor; four_velocity=("vel1","vel2"))
            @test_throws ArgumentError quantity(s,model,:lorentz_factor; four_velocity=("vel1","vel2","vel1"))
            @test_throws ArgumentError quantity(s,model,:comoving_magnetic_pressure; magnetic=("missing","Bcc2","Bcc3"))
            @test_throws ArgumentError quantity(s,model,:lorentz_factor; stored="missing")
            @test_throws ArgumentError quantity(s,model,:comoving_magnetic_pressure;
                                                   magnetic=("Bcc1","Bcc2","Bcc3"), stored="pB_fluid")
            @test_throws ArgumentError quantity(s,model,:lorentz_factor; magnetic=("Bcc1","Bcc2","Bcc3"))
            @test_throws DomainError read_field(s,quantity(s,model,:lorentz_factor; stored="press"),1)
            @test_throws DomainError read_field(s,quantity(s,model,:lab_velocity_2; stored="rho"),1)
        end
        h5open(path,"r+") do h
            h["extras"][1,1,1,1,2] = NaN  # vel1
            h["extras"][1,1,1,2,5] = Inf  # Bcc1
            h["extras"][1,1,1,1,8] = -1.  # invalid stored pressure
        end
        open_snapshot(path) do s
            @test_throws DomainError read_field(s,quantity(s,model,:lorentz_factor),1)
            @test_throws DomainError read_field(s,quantity(s,model,:comoving_magnetic_pressure),2)
            q = quantity(s,model,:comoving_magnetic_pressure; stored="pB_fluid")
            @test_throws DomainError read_field(s,q,1)
            @test first(read_field(s,q,2)) == pB
            @test first(read_field(s,quantity(s,model,:lorentz_factor),2)) == 1
        end
        h5open(path,"r+") do h
            HDF5.delete_attribute(h,"Coordinates")
            attributes(h)["Coordinates"] = "kerr-schild"
        end
        open_snapshot(path) do s
            @test_throws ArgumentError quantity(s,model,:lorentz_factor)
        end
    end
end
