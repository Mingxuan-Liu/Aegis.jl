"""A source field paired with optional, explicitly user-supplied labels."""
struct FieldInterpretation
    source::FieldRef
    alias::Union{Nothing,String}
    units::Union{Nothing,String}
    description::Union{Nothing,String}
end

"""
    interpret_field(snapshot, selector; alias=nothing, units=nothing, description=nothing)

Attach user-supplied labels to a field from this snapshot. The returned descriptor
is also a selector for `read_field` and `summarize`. Keep several descriptors in
a named tuple or dictionary to look them up by your own aliases. `source`
preserves the original dataset, slot and name. Unspecified labels stay `nothing`.
Labels describe the stored numbers; no unit conversion, rescaling or change of
centering is performed.
"""
function interpret_field(snapshot::Snapshot, selector; alias=nothing, units=nothing, description=nothing)
    return FieldInterpretation(resolve_field(snapshot, selector), alias, units, description)
end

resolve_field(snapshot::Snapshot, field::FieldInterpretation) = resolve_field(snapshot, field.source)
