"""
Native MeshBlock analysis of Athena++ fluid output, including Newtonian and
special-relativistic primitives, Cartesian derivatives and Fourier spectra.
See the package README for a verified workflow and the supported data contract.
"""
module Aegis

using HDF5

export open_snapshot, discover_snapshots, resolve_field, read_field, block_coordinates
export interpret_field
export read_athinput, input_blocks, parameter_names, parameter
export read_history, history_column, restart_rows
export CartesianMesh, Region, selectedblocks, cell_volume, summarize
export NewtonianMHD, IdealGas, Isothermal, quantity
export SpecialRelativisticMHD, magnetic_fluctuation
export BinAxis, cell_map, distribution, conditional_stats, slice_plane, project, export_product
export DifferentialGrid, gradient, vector_gradient, divergence, curl
export shell_spectrum, spectrum_workspace

include("athdf.jl")
include("interpretation.jl")
include("textreaders.jl")
include("reductions.jl")
include("quantities.jl")
include("relativistic.jl")
include("differentials.jl")
include("products.jl")
include("fluctuations.jl")
include("spectra.jl")

end
