# The three private producer configurations were inspected individually. The
# oracle samples all raw components via h5py in C storage order.
function check_private_quantities(s, input, ref, closure)
    fields = Dict(zip(ref["names"],ref["samples"]))
    r = fields["rho"]
    v2 = sum(fields["vel$i"]^2 for i in 1:3)
    b2 = sum(fields["Bcc$i"]^2 for i in 1:3)
    if closure === :isothermal
        cs = parameter(input,"hydro","iso_sound_speed",Float64)
        p = r*cs^2
        model = NewtonianMHD(eos=Isothermal(cs))
    else
        gamma = parameter(input,"hydro","gamma",Float64)
        p = fields["press"]
        cs = sqrt(gamma*p/r)
        model = NewtonianMHD(eos=IdealGas(gamma))
    end
    expected = Dict(:speed=>sqrt(v2), :magnetic_magnitude=>sqrt(b2),
        :kinetic_energy_density=>r*v2/2, :magnetic_energy_density=>b2/2,
        :magnetic_pressure=>b2/2, :gas_pressure=>p, :plasma_beta=>2p/b2,
        :sound_speed=>cs, :alfven_speed=>sqrt(b2/r),
        :sonic_mach=>sqrt(v2)/cs, :alfven_mach=>sqrt(r*v2/b2))
    if closure === :isothermal
        @test_throws ArgumentError quantity(s,model,:thermal_energy_density)
    else
        expected[:thermal_energy_density] = p/(model.eos.gamma-1)
    end
    @test_throws ArgumentError quantity(s,model,:temperature)
    ranges = Tuple(i:i for i in ref["indices"])
    for (name,want) in expected
        value = only(read_field(s,quantity(s,model,name),ref["block"]; ranges))
        @test isapprox(value,want; rtol=8eps(Float64))
    end
end
