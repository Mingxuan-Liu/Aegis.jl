"""
Native MeshBlock access and Cartesian reductions for Athena++ mesh output.
See the package README for a verified workflow and the supported data contract.
"""
module Aegis

using HDF5

export open_snapshot, discover_snapshots, resolve_field, read_field, block_coordinates
export interpret_field
export read_athinput, input_blocks, parameter_names, parameter
export read_history, history_column, restart_rows
export CartesianMesh, Region, selectedblocks, cell_volume, summarize

include("athdf.jl")
include("interpretation.jl")
include("textreaders.jl")
include("reductions.jl")

end
