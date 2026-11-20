# Reuse the independent h5py sample in integration.jl. The specified axis
# checks the calculation on each output and makes no claim about its initial field.
function check_private_fluctuations(s, ref)
    values = Dict(zip(ref["names"], ref["samples"]))
    bx, by, bz = (values["Bcc$i"] for i in 1:3)
    reference = (0.25,-0.125,0.0625)
    # Resolve along (3/5,4/5,0), and independently onto the two
    # transverse axes (-4/5,3/5,0) and (0,0,1).
    x,y,z = (bx,by,bz) .- reference
    expected = ((3x/5+4y/5)^2/2, ((3y/5-4x/5)^2+z^2)/2)
    ranges = Tuple(i:i for i in ref["indices"])
    for (component, answer) in zip((:parallel,:transverse), expected)
        q = magnetic_fluctuation(s,NewtonianMHD(),component; reference,direction=(3,4,0))
        value = only(read_field(s,q,ref["block"]; ranges))
        @test isapprox(value,answer; rtol=16eps(), atol=4eps()*(x^2+y^2+z^2))
    end
end
