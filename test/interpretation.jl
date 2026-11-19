@testset "User descriptions preserve file identity and numbers" begin
    mktempdir() do dir
        files = example_run(dir)
        open_snapshot(files.snapshot) do s
            aliases = (sample=interpret_field(s, "tracer"; alias="sample", units="arb."),)
            f = aliases.sample
            @test f.source.name == "tracer"
            @test f.alias == "sample"
            @test f.units == "arb."
            @test isnothing(f.description)
            @test read_field(s, f, 1) == read_field(s, "tracer", 1)
            result = summarize(s, f)
            @test result.integral == 819
            @test result.provenance.field_selector === f
            @test result.provenance.field == f.source
        end
    end
end
