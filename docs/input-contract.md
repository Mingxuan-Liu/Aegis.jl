# Native Athena++ MHD mesh output

## Evidence and units

A snapshot records the mesh and sampled fields at its header time. The column names and printed times in a history belong to its own output stream. An `athinput` describes a requested run configuration. It can predate command-line overrides, a restart, or a different executable. Read these sources separately. Use the snapshot's stored arrays, field list, geometry and time for work on that snapshot. Resolve any conflicting interpretation against the producer before applying a physical conversion.

For the nonrelativistic Athena++ MHD equations, magnetic permeability is one. The magnetic pressure and energy density are both $B^2/2$ in code units. The ideal-gas MHD energy density, when the adiabatic model applies, includes internal, bulk kinetic and magnetic energy. The isothermal model has a different closure. The existence of a `gamma` input entry alone cannot identify either the compiled EOS or a stored pressure field. Density, velocity and length scales for conversion to physical units depend on the run. For example, in Gaussian CGS one would use a magnetic scale $B_0=\sqrt{4\pi\rho_0}\,v_0$ given independently established mass-density and velocity scales. The package returns values in the stored units. A user-supplied unit label records the user's interpretation. See Stone et al. (2020), §§2.1 and 3.1, especially Eqs. (6)–(8). [Paper](https://doi.org/10.3847/1538-4365/ab929b)

The MHD-PIC equations evolve the thermal fluid and energetic particles separately. With backreaction, the gas momentum and energy equations contain particle exchange terms. A gas energy field cannot supply the combined gas–particle energy. The artificial particle light speed controls the momentum–velocity and kinetic-energy relations of the particles. Applying it as a fluid velocity conversion, or importing a different reduced-speed-of-light prescription, requires additional physical assumptions. An expanding-box calculation also needs metric information before computing physical volumes. See Sun & Bai (2023), §§2.1–2.3, Eqs. (2)–(6), footnote 1, and §3 (expanding box). [Paper](https://arxiv.org/abs/2304.10568)

## From file axes to a physical cell

| Evidence in the file | Interpretation for supported outputs |
| --- | --- |
| `DatasetNames`, `NumVariables`, flattened `VariableNames` | The counts partition the names in order among the listed datasets. A dataset and its local variable slot identify a source, even if a label occurs twice. An uninterpreted custom field can still be read. |
| Five-dimensional field dataset | HDF5/C order is (variable, block, x3, x2, x1). With HDF5.jl, a selected field and block yields a Julia array indexed `[i,j,k]` in physical (x1,x2,x3) order, including a length-one axis in 2D. File values retain the precision chosen for output. |
| `Levels`, `LogicalLocations` | Each stored block has a refinement level relative to the root mesh, and three zero-based logical indices at that level. Julia block indices identify storage slots. Neighboring slots need not be spatial neighbors. |
| `x1f` … `x3f`, `x1v` … `x3v` | Per-block faces and volume-center coordinates locate the cells. For Cartesian cells, multiply the three successive face differences, including the recorded width in a singleton direction. This gives a volume in code length cubed. |
| `RootGridSize`, `MeshBlockSize`, `RootGridX1` … `RootGridX3`, `Coordinates`, `MaxLevel` | Root cell counts, output extents per block, bounds and mesh-spacing parameters constrain coverage. Some Athena++ writes include ghosts or collapse axes for sums or slices. Compare the stored mesh against a complete nonoverlapping native tiling before interpreting a full domain. |
| `Time`, `NumCycles` | These describe the state at output. An output number indexes a stream. Matching numbers alone cannot synchronize two streams, resets, or histories printed at a different precision. |

**Centering.** The `prim`/`cons` MHD output path in the inspected writers takes `Bcc1`, `Bcc2`, `Bcc3` from the solver's cell-centered magnetic array. Explicit `b1`, `b2`, `b3` output takes values from the corresponding lower faces. The writers can place *both* in arrays with one entry per cell and even write `Center="Cell"` in XDMF. There is no upper boundary face in the packed face-component output. A volume statistic of a cell-centered field has a different meaning from a statistic of those face samples. Aegis treats the labels `B1`, `B2` and `B3` as ambiguous in this format. Raw reads can retrieve them. A custom output's stored sample and the physical quantity assigned to it can likewise differ. Check the output selection and producer for custom conventions.

**Coverage.** Geometry-aware operations target complete, unrefined Cartesian native cells. They use the saved faces, including nonuniform spacing. A center-based region includes each selected cell's full represented volume. An exact truncated-cell integral requires a reconstruction or intersection rule. Missing blocks, refinement, ghost overlap, sums, slices of higher-dimensional grids and alternate coordinates need their own coverage/metric rules. Raw block reads preserve their native organization when the storage layout is compatible. A field's reference frame, EOS, dimensional conversion, pgen formula or physical density role must be established by the caller when needed for interpretation. Native field reads return the stored components in their original frame. The [SR diagnostics](relativity.md) specify which lab-frame inputs determine each rest-frame quantity.

## Supplied integration cases

The following inventory describes the local examples. Field labels and run identifiers in this table are validation evidence, never reader requirements. The inspected snapshot headers contain no producing executable commit. The source revisions below establish compatibility, while the build identity of the files remains unknown.

| Local run | Snapshots in stream, code time | Root cells, cells per block, number of blocks | Stored groups and fields | History and input |
| --- | --- | --- | --- | --- |
| 3D MHD, 128 | 11, 0 to 10 | 128³, 32³, 64 | `hydro` (16) | 1,001 rows × 26 columns, one input |
| 3D MHD, 256 | 1, 10 | 256³, 32³, 512 | `hydro` (16) | 1,001 rows × 26 columns, one input |
| 2D MHD-PIC | 13, 0 to 600 | 1024 × 4096 × 1, 64 × 64 × 1, 1,024 | `prim` (4) and `B` (3) | 1,201 rows × 50 columns, one input |

All three report `cartesian`, level zero and Float32 field and coordinate data. The 3D files contain `rho`, `press`, three velocities, three cell-centered magnetic components, and eight user fields (`force_x`, `force_y`, `force_z`, `cooling`, `current_density`, `ohmic_heating`, `viscous_heating`, `conductive_heating`). Both 3D domains extend from −4 to 4 along each axis. The final snapshots record cycles 10,167 and 25,456, respectively. Both inputs request `prim,uov`. Their histories contain names with spaces. The 128³ run's job log supplies run-specific dimensional scales. Reuse on another run requires independent verification.

The 2D input records an isothermal MHD-PIC configuration and requests `prim`. The output contains fluid density, three velocities, and three cell-centered magnetic components. Pressure and particle snapshots are absent. The history includes user columns and some adjacent indexed labels with no space between them. At the final time the snapshot records cycle 1,204,717. The x1 and x2 faces span approximately ±15.707963 and ±62.831852, respectively. The full x3 interval is [−0.5, 0.5] in code units. Multiplication by this thickness gives the volume of the represented slab.

**Independent mass check.** The standard history writer accumulates conserved gas density times active-cell volume. For these nonrelativistic Cartesian cases, the last history row and last supplied mesh file have the same printed time. Independent h5py reads of the saved density, integrated one block at a time with Float64 differences of saved faces, gave the following totals in code mass units. Rounding the totals to the history's six significant digits reproduces its printed mass. This single-state comparison places no constraint on evolution between snapshots.

| Run | Time | Integrated mesh density | Solver history mass as printed |
| --- | ---: | ---: | ---: |
| 128³ MHD | 10 | 511.9999999748 | 5.12000e+02 |
| 256³ MHD | 10 | 512.0000000020 | 5.12000e+02 |
| 2D MHD-PIC | 600 | 3947.841619869 | 3.94784e+03 |

## Implementation references

- The compatible MHD-PIC source shows variable assignment in `OutputType::LoadOutputData`, extent/axis serialization in `ATHDF5Output::WriteOutputFile`, gas mass in `HistoryOutput::WriteOutputFile`, and Cartesian volume in `Coordinates::CellVolume`. These symbols distinguish solver staggering from stored array dimensions. [Field assignments](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/outputs.cpp) · [HDF5 writer](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/athena_hdf5.cpp) · [History](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/history.cpp) · [Volumes](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/coordinates/coordinates.cpp)
- The MHD fork supports combined output selections such as `prim,uov`. Its parser implements block/key source syntax independently of the scientific meaning of the value. Its reference reader's `athdf(..., raw=True)` path retains native block arrays without regridding. [Output selection](https://github.com/dfielding14/athena-1/blob/d2aa0b0ffc3b6b215721693cf82946cfa5582616/src/outputs/outputs.cpp) · [Parameter input](https://github.com/dfielding14/athena-1/blob/d2aa0b0ffc3b6b215721693cf82946cfa5582616/src/parameter_input.cpp) · [Reference reader](https://github.com/dfielding14/athena-1/blob/d2aa0b0ffc3b6b215721693cf82946cfa5582616/vis/python/athena_read.py)
