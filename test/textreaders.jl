using Test

@testset "Athena++ input syntax and typed access" begin
    # Repeated block, repeated key, tabs, arbitrary text, empty fields, end marker.
    data = "  # comments and blank lines\r\n\n<custom>\n" *
           " mode\t=  a=b, c  # source comment\n" *
           "rate = 4E-3\n<empty>\n<custom> # repeat\n" *
           "rate = +2.5e-2  # revised\nflag = tRuE\nzero = 0\n" *
           "one = 1\noff = FALSE\nnum = -02\nempty = # intentionally blank\n" *
           "expr = sqrt(4)\nliteral = a&\n<par_end>\ninvalid trailer\n"
    input = read_athinput(IOBuffer(data); source="manufactured.input")
    @test input_blocks(input) == ["custom", "empty"]
    @test parameter_names(input, "empty") == []
    @test parameter_names(input, "custom") == ["mode","rate","flag","zero","one","off","num","empty","expr","literal"]
    @test parameter(input, "custom", "mode") == "a=b, c"
    @test parameter(input, "custom", "rate") == "+2.5e-2"
    @test parameter(input, "custom", "rate", Float64) == 0.025
    @test parameter(input, "custom", "num", Int) == -2
    @test parameter(input, "custom", "flag", Bool)
    @test !parameter(input, "custom", "zero", Bool)
    @test parameter(input, "custom", "one", Bool)
    @test !parameter(input, "custom", "off", Bool)
    @test parameter(input, "custom", "empty") == ""
    @test parameter(input, "custom", "expr") == "sqrt(4)"
    @test parameter(input, "custom", "literal") == "a&"
    @test [e.line for e in input.entries if e.name == "rate"] == [5,8]
    @test input.entries[1].text == " mode\t=  a=b, c  # source comment"
    @test_throws ArgumentError parameter(input, "custom", "expr", Float64)
    @test_throws ArgumentError parameter(input, "custom", "rate", Int)
    @test_throws ArgumentError parameter(input, "custom", "num", Bool)
    @test_throws ArgumentError parameter(input, "CUSTOM", "flag")
    @test_throws ArgumentError parameter_names(input, "absent")
    for invalid in ["x=1", "<unfinished\n", "  <b>\nx=2", "<b>\n=4", "<b>\nx=as &\ncontinued text"]
        @test_throws ArgumentError read_athinput(IOBuffer(invalid))
    end
    # The verified branches remove tabs, even inside lexical tokens.
    @test parameter(read_athinput(IOBuffer("<b>\nk=e\tf\n")), "b", "k") == "ef"
end

@testset "History schema identity, labels and branch selection" begin
    data = """
    # Athena++ history data
    # [1]=time [2]=Mean speed[3]=flux α [4]=flux α
    0 10 20 30
    1 11 21 31
    3 13 23 33
    # Athena++ history data
    # [1]=time [2]=Mean speed[3]=flux α [4]=flux α
    1 12 22 32
    2 14 24 34
    # Athena++ history data
    # [1]=added[2]=Mean speed[3]=time
    -1 15 2
    -2 16 4
    """
    h = read_history(IOBuffer(data); source="manufactured.hst")
    @test h.schemas[1].names == ["time", "Mean speed", "flux α", "flux α"]
    @test h.schemas[1].names == h.schemas[2].names
    @test h.schemas[3].names == ["added", "Mean speed", "time"]
    @test h.schemas[2].line == 7
    @test [r.schema for r in h.rows] == [1,1,1,2,2,3,3]
    @test h.rows[4].line == 8 && h.rows[4].text == "1 12 22 32"
    @test history_column(h, "time") == [0,1,3,1,2,2,4]
    @test restart_rows(h) == [1,4,6,7]
    @test history_column(h, "Mean speed"; rows=restart_rows(h)) == [10,12,15,16]
    @test history_column(h, 4; schema=1) == [30,31,33]
    @test history_column(h, "added"; schema=3, rows=restart_rows(h)) == [-1,-2]
    @test_throws ArgumentError history_column(h, 1)
    @test_throws ArgumentError history_column(h, 5; schema=1)
    @test_throws ArgumentError history_column(h, "flux α"; schema=1)
    @test_throws ArgumentError history_column(h, "added")
    @test_throws ArgumentError history_column(h, "time"; schema=4)
    # Missing time or nonfinite samples remain available until time resolution.
    unlabelled = read_history(IOBuffer("# [1]=[2]=v [3]=v\n0 NaN Inf\n"))
    @test isempty(unlabelled.schemas[1].names[1])
    @test isinf(only(history_column(unlabelled, 3)))
    @test_throws ArgumentError restart_rows(unlabelled)
    nonfinite = read_history(IOBuffer("# [1]=time\nNaN\n"))
    @test isnan(only(history_column(nonfinite, "time")))
    @test_throws ArgumentError restart_rows(nonfinite)
    duplicate_times = read_history(IOBuffer("# [1]=time [2]=time\n0 1\n"))
    @test_throws ArgumentError restart_rows(duplicate_times)
    empty_table = read_history(IOBuffer("# [1]=time\n"))
    @test isempty(restart_rows(empty_table))
    @test isempty(history_column(empty_table, "time"; schema=1))
    @test_throws ArgumentError history_column(empty_table, "unknown"; schema=1)
    @test_throws ArgumentError history_column(empty_table, "unknown")
    for invalid in ["", "1 2", "# [2]=v\n0", "# [1]=a [3]=b\n0 1",
                    "# [1]=a [2]=b\n0", "# [1]=a\ninvalid",
                    "# [1]=a\n0\n# Athena++ history data\n1",
                    "# Athena++ history data\n"]
        @test_throws ArgumentError read_history(IOBuffer(invalid))
    end
end
