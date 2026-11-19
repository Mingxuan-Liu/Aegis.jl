# Input and history grammars verified in PrincetonUniversity/MHD-PIC-in-Athena-
# 3eb8ec0eeba5e213a3c0f23133fe6f1f287e2513, src/parameter_input.cpp
# (LoadFromStream, ParseLine, AddParameter) and src/outputs/history.cpp
# (HistoryOutput::WriteOutputFile). The inspected MHD writer also uses indexed headers.

struct InputEntry
    block::String
    name::String
    value::String
    line::Int
    text::String
end

struct AthenaInput
    source::String
    blocks::Vector{String}
    entries::Vector{InputEntry}
    latest::Dict{Tuple{String,String},Int}
end

_text_error(source, line, message) = throw(ArgumentError("$source:$line: $message"))

"""
    read_athinput(path)
    read_athinput(io; source="<io>")

Read case-sensitive Athena++ `<block>` and `key = value # comment` syntax.
Repeated blocks refer to the same block. Each assignment in `input.entries`
has its original line number, text and parsed value. Lookup uses the last
assignment. No defaults, environment expansion, units or EOS are inferred.

The parser stops at `<par_end>` or EOF. In the supported Athena++ and MHD-PIC
dialects, tabs are discarded before parsing, a hash starts a comment even
inside quotes, and each assignment occupies one line. Backslashes and
ampersands in values are literal. Block headers must start in column one
after tab removal: the verified C++ parser includes the closing bracket in
the stored name when a header is indented with spaces. See
[ParameterInput::LoadFromStream and ParseLine](https://github.com/PrincetonUniversity/MHD-PIC-in-Athena-/blob/3eb8ec0eeba5e213a3c0f23133fe6f1f287e2513/src/parameter_input.cpp).
"""
read_athinput(path::AbstractString) = open(io -> read_athinput(io; source=abspath(path)), path)

function read_athinput(io::IO; source::AbstractString="<io>")
    blocks = String[]
    entries = InputEntry[]
    latest = Dict{Tuple{String,String},Int}()
    block = nothing
    for (line_no, raw) in enumerate(eachline(io))
        untabbed = replace(raw, '\t' => "")
        line = strip(first(split(untabbed, '#'; limit=2)))
        isempty(line) && continue
        startswith(line, "<par_end>") && break
        if startswith(line, "<")
            startswith(untabbed, "<") ||
                _text_error(source, line_no, "left-align the <block> header; Athena++ misreads indented block names")
            m = match(r"^<([^<>]+)>$", line)
            isnothing(m) && _text_error(source, line_no, "expected a <block> header")
            block = String(m.captures[1])
            block in blocks || push!(blocks, block)
            continue
        end
        isnothing(block) && _text_error(source, line_no, "assignment before the first <block>")
        pair = split(line, '='; limit=2)
        length(pair) == 2 || _text_error(source, line_no, "expected key = value in <$block>")
        key, value = strip.(pair)
        isempty(key) && _text_error(source, line_no, "missing parameter name in <$block>")
        push!(entries, InputEntry(block, key, value, line_no, raw))
        latest[(block, key)] = length(entries)
    end
    return AthenaInput(String(source), blocks, entries, latest)
end

"""Return block names in order of their first appearance, including empty blocks."""
input_blocks(input::AthenaInput) = copy(input.blocks)

"""Return arbitrary parameter names in a block in order of their first assignment."""
function parameter_names(input::AthenaInput, block::AbstractString)
    block in input.blocks || throw(ArgumentError("$(input.source): no <$block> block"))
    return unique([e.name for e in input.entries if e.block == block])
end

function _input_entry(input, block, name)
    index = get(input.latest, (block, name), 0)
    index > 0 || throw(ArgumentError("$(input.source): no parameter <$block>/$name"))
    return input.entries[index]
end

"""
    parameter(input, block, name[, T])

Return the last assigned string or convert the entire value to a decimal
integer, floating-point number or Boolean. Boolean values are case-insensitive
`true`/`false` or `1`/`0`. Missing keys and failed conversions raise an error
with source information. Strict conversion prevents malformed inputs becoming
plausible numbers through the solver's permissive C numeric prefix conversion.
For a custom interpretation, retrieve the string and apply a downstream
function. No physical scale or role is attached to any key.
"""
parameter(input::AthenaInput, block::AbstractString, name::AbstractString, ::Type{String}=String) =
    _input_entry(input, block, name).value

function parameter(input::AthenaInput, block::AbstractString, name::AbstractString,
                   ::Type{T}) where {T<:Union{Integer,AbstractFloat}}
    entry = _input_entry(input, block, name)
    value = entry.value
    parsed = if T === Bool
        token = lowercase(value)
        token in ("true", "1") ? true : token in ("false", "0") ? false : nothing
    elseif T <: Integer
        tryparse(T, value; base=10)
    else
        tryparse(T, value)
    end
    isnothing(parsed) && _text_error(input.source, entry.line,
        "cannot convert <$block>/$name = $(repr(value)) to $T")
    return parsed
end

struct HistorySchema
    line::Int
    names::Vector{String}
    text::String
end

struct HistoryRow
    line::Int
    schema::Int
    values::Vector{Float64}
    text::String
end

struct AthenaHistory
    source::String
    schemas::Vector{HistorySchema}
    rows::Vector{HistoryRow}
end

# Athena writes "[i]=label" with padding but no required inter-column delimiter.
# Match positions delimit UTF-8 names, including embedded spaces and punctuation.
function _history_schema(line, raw, source, line_no)
    markers = collect(eachmatch(r"\[(\d+)\]\s*=\s*", line))
    isempty(markers) && _text_error(source, line_no, "expected indexed history columns [1]=...")
    prefix = strip(SubString(line, firstindex(line), prevind(line, first(markers).offset)))
    isempty(prefix) || _text_error(source, line_no, "unexpected text before first column")
    names = String[]
    for (i, marker) in enumerate(markers)
        tryparse(Int, marker.captures[1]) == i ||
            _text_error(source, line_no, "expected consecutive column numbers starting at 1; expected [$i]")
        start = marker.offset + ncodeunits(marker.match)
        stop = i == length(markers) ? lastindex(line) : prevind(line, markers[i+1].offset)
        push!(names, strip(SubString(line, start, stop)))
    end
    return HistorySchema(line_no, names, raw)
end

"""
    read_history(path)
    read_history(io; source="<io>")

Read every numbered `.hst` header and data row in source order. Each member of
`history.schemas` records the header line, its original text and `names` in
positional order, including duplicates. Each member of `history.rows` records
the source line, raw text, schema index and numeric `values::Vector{Float64}`.
Column names may contain spaces or touch the next numbered marker. Numbers
are parsed at the precision printed by the writer. Units, weights and
definitions of custom columns must come from the user or the simulation.

Repeated headers each start a schema record, even for identical columns.
Repeated and decreasing times, missing time labels and nonfinite data remain
visible. Malformed numeric records and column counts raise errors. This reader
holds the text history in memory, with space proportional to the input size.
The indexed header convention is from
[HistoryOutput::WriteOutputFile](https://github.com/PrincetonUniversity/MHD-PIC-in-Athena-/blob/3eb8ec0eeba5e213a3c0f23133fe6f1f287e2513/src/outputs/history.cpp).
"""
read_history(path::AbstractString) = open(io -> read_history(io; source=abspath(path)), path)

function read_history(io::IO; source::AbstractString="<io>")
    schemas = HistorySchema[]
    rows = HistoryRow[]
    pending_header = false
    for (line_no, raw) in enumerate(eachline(io))
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#")
            comment = strip(chopprefix(line, "#"))
            if comment == "Athena++ history data"
                pending_header && _text_error(source, line_no, "previous descriptor has no column header")
                pending_header = true
            elseif startswith(comment, "[")
                push!(schemas, _history_schema(comment, raw, source, line_no))
                pending_header = false
            end
            continue
        end
        (isempty(schemas) || pending_header) &&
            _text_error(source, line_no, "history row has no indexed column header")
        tokens = split(line)
        expected = length(last(schemas).names)
        length(tokens) == expected ||
            _text_error(source, line_no, "expected $expected values, found $(length(tokens))")
        values = Vector{Float64}(undef, expected)
        for (i, token) in enumerate(tokens)
            value = tryparse(Float64, token)
            isnothing(value) && _text_error(source, line_no, "invalid number in column $i: $token")
            values[i] = value
        end
        push!(rows, HistoryRow(line_no, length(schemas), values, raw))
    end
    pending_header && throw(ArgumentError("$source: final descriptor has no column header"))
    isempty(schemas) && throw(ArgumentError("$source: no indexed history column header"))
    return AthenaHistory(String(source), schemas, rows)
end

function _column_index(source, schema, name::AbstractString)
    idx = findfirst(==(name), schema.names)
    isnothing(idx) && _text_error(source, schema.line, "missing history column $(repr(name))")
    !isnothing(findnext(==(name), schema.names, idx + 1)) &&
        _text_error(source, schema.line, "multiple columns named $(repr(name)); select by index and schema")
    return idx
end

function _column_index(source, schema, index::Integer)
    checkbounds(Bool, schema.names, index) ||
        _text_error(source, schema.line, "column position $index is out of bounds")
    return Int(index)
end

"""
    history_column(history, selector; schema=nothing, rows=eachindex(history.rows))

Return Float64 values for the named column across the selected source rows.
Lookup follows the exact name even if a later header changes its position.
Names must resolve uniquely in every included header, even an empty one.
For positional lookup
across different layouts, supply the 1-based `schema` keyword. Specifying a
schema restricts column lookup and the result to that header. `rows` filters
the values after column lookup and accepts an iterable of raw row indices,
including the output of [`restart_rows`](@ref).
"""
function history_column(history::AthenaHistory, selector::Union{AbstractString,Integer};
                        schema::Union{Nothing,Integer}=nothing, rows=eachindex(history.rows))
    if !isnothing(schema) && !checkbounds(Bool, history.schemas, schema)
        throw(ArgumentError("$(history.source): unknown schema $schema"))
    end
    eligible = isnothing(schema) ? eachindex(history.schemas) : (schema,)
    positions = zeros(Int, length(history.schemas))
    for s in eligible
        header = history.schemas[s]
        if selector isa Integer && isnothing(schema) && header.names != first(history.schemas).names
            _text_error(history.source, header.line, "positional selection across layouts requires schema")
        end
        positions[s] = _column_index(history.source, header, selector)
    end
    values = Float64[]
    for i in rows
        row = history.rows[i]
        index = positions[row.schema]
        iszero(index) && continue
        push!(values, row.values[index])
    end
    return values
end

"""
    restart_rows(history; time="time")

Return raw row indices of the final increasing time branch. On each row at
time `t`, discard previously retained rows with time ≥ `t`, then append the
new row. This opt-in policy treats a time reversal or tie (including ties
caused by printed precision) as a branch replacement. No rows are sorted
and the original history stays available. Header changes alone do not
truncate a branch. Each represented schema must have one exact, finite time
column under the supplied name; the numeric times must share a scale. The
algorithm uses linear time and O(number of rows + schemas) extra space.
"""
function restart_rows(history::AthenaHistory; time::AbstractString="time")
    retained = Int[]
    times = Float64[]
    positions = Dict{Int,Int}()
    for (i, row) in enumerate(history.rows)
        index = get!(positions, row.schema) do
            _column_index(history.source, history.schemas[row.schema], time)
        end
        t = row.values[index]
        isfinite(t) || _text_error(history.source, row.line, "nonfinite restart time")
        while !isempty(times) && last(times) >= t
            pop!(times)
            pop!(retained)
        end
        push!(times, t)
        push!(retained, i)
    end
    return retained
end
