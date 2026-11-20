# Native distributions and spatial products

These operations use the same verified Cartesian mesh and half-open center selection as `summarize`. They accept stored fields, user-labelled fields, model-specific fluid quantities, guide-field fluctuation energies, and scalar expressions defined with `cell_map`. A region contributes whole native cells with their stored face widths. Code lengths, volumes and fields refer to the snapshot frame and the units described in the [input contract](input-contract.md).

## Binning with a specified measure

Supply bin edges in the physical coordinate of each scalar, even for a per-decade distribution. By default each cell contributes $W_i=\Delta V_i$. A weight selector gives $W_i=w_i\Delta V_i$, for example a density selector gives mass. The special weight `:count` gives $W_i=1$. Weights must be finite and nonnegative. A histogram over a selected domain exposes the following quantities.

| Quantity | Definition |
| --- | --- |
| `totals[J]` | Sum of cell weights in bin $J$, with joint axes in argument order |
| `included_weight` | Sum of the bin totals |
| `excluded_weight` | Sum of selected weights outside one or more axes, counted once per cell |
| `widths[d][j]` | $e_{j+1}-e_j$ for `:linear`, or $\log_{10}(e_{j+1}/e_j)$ decades for `:dex` |
| `pdf[J]` | Bin total divided by the included weight and the product of the axis widths |

All bins include their left edge. The final bin on each axis includes its right edge. An observation exactly on an internal edge enters the bin to its right. A positive included weight gives $\sum_J p_J\prod_d \Delta u_{d,J}=1$. Empty bins then have zero density. With zero included weight, every density is NaN and the bin totals and excluded weight remain available. Nonpositive finite coordinates on a dex axis fall outside its positive edges. Nonfinite coordinates and weights raise errors. The edge and density conventions follow the NumPy histogram definition, applied independently per axis. [NumPy histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)

A conditional statistic returns $N_j=\sum_i W_i q_i$, $D_j=\sum_i W_i$, and $N_j/D_j$ within each bin. Empty and zero-weight bins have NaN means. For heating and cooling fields $H,C$, the conditional mean of the local ratio uses `cell_map(/, s, H, C; name="H/C", units="1")` as the response, with a domain in which $C$ is finite and nonzero. The ratio of integrated powers divides the separate conditional numerators for $H$ and $C$ with volume weighting. The two definitions differ when the local cooling power varies.

## Plane and line-of-sight geometry

`slice_plane` takes the cell containing the coordinate `at`. At a shared interface it takes the layer on the positive side. The domain's upper face lies outside this half-open convention. If the chosen layer's center lies outside the selected region, the slice raises an error. Each pixel covers exactly one native transverse cell and retains the field's units.

`project` follows the ray along `axis`. A line integral gives $P=\sum_i q_i\Delta l_i$. The path mean divides $P$ by the represented path length. A weighted mean gives $\sum_i w_iq_i\Delta l_i / \sum_i w_i\Delta l_i$ for each ray. A zero denominator gives NaN. The transverse axes are the remaining physical dimensions in increasing order. The two face vectors give pixel areas $A_{ab}=\Delta u_a\Delta v_b$, so $\sum_{ab}A_{ab}P_{ab}$ recovers the volume integral for the same selected field and region. This remains true with unequal native widths and a physical singleton thickness. [yt projection definitions](https://yt-project.org/doc/visualizing/plots.html)

For the tested two-layer slab, the path lengths are 1 and 3, $q=(2,6)$ and $\rho=(3,1)$. The column integral is $2+18=20$, its path mean is 5 and its mass-weighted mean is $(6+18)/(3+3)=4$. The native test splits the same slab across shuffled MeshBlocks and compares all three normalizations.

## Runnable example

Run the following from the checkout root. The fixture writes into a temporary directory.

```julia
using Aegis
include("test/fixtures.jl")
mktempdir() do dir
    run = example_run(dir)
    open_snapshot(run.snapshot) do s
        mesh = CartesianMesh(s)
        tracer = interpret_field(s, "tracer"; units="tracer unit")
        rho = interpret_field(s, "rho"; units="code density")
        x = BinAxis(tracer, [1, 2, 4, 6]; measure=:dex)
        y = BinAxis(rho, [2, 4, 7])
        joint = distribution(mesh, x, y; weight=rho)
        @assert joint.included_weight == 1008
        means = conditional_stats(mesh, x, rho)
        plane = slice_plane(mesh, tracer; axis=2, at=1)
        column = project(mesh, tracer; axis=1)
        area = diff(column.faces[1]) .* diff(column.faces[2])'
        @assert sum(column.values .* area) == 819
        export_product(joinpath(dir,"joint.h5"), joint)
        export_product(joinpath(dir,"column.h5"), column)
    end
end
```

## Check on a simulation output

The final 128³ filament snapshot (code time 10) contains the stored field `ohmic_heating`. For the region $[0,1)\times[-2,1)\times[-1,2)$, a mass histogram with physical edges $0$, $0.01$, $0.1$, and $10^{20}$ gives bin masses 0.0309128562, 0.8359968167 and 10.4263247542 in code mass units. Integrating the field along $x_2$ and summing the map times its pixel areas gives 3.0088604754 in field-times-volume units. The $x_2=0.5$ slice has a pixel sum of 374.2504986. Independent NumPy calculations from the stored HDF5 slabs, faces, and centers reproduced these checks. The example attests to access and weighting of the pgen field; its stored name alone does not define a dimensional heating conversion.

## Output format and resource use

`export_product` writes a new HDF5 file. Properties of the returned Julia named tuple become root datasets or named groups. Tuples become groups with one-based numbered children. Absent optional properties are omitted. Root attributes record `aegis_product_schema` (currently 1) and the package version. The `metadata` group records the source path, snapshot time, available cycle, region, coordinate and weighting conventions, original field identities, and known units. Mapped fields carry their declared labels and input descriptors. Callers should retain their source code for custom functions. HDF5.jl preserves Julia physical-axis ordering on readback. In an HDF5/C or h5py reader, array dimensions appear in reverse order. For a 2D product, transposing the array restores the axis order given in the numbered coordinate groups.

Histograms, conditional statistics and projections stream one block slab at a time. Top-level repeated field selectors share that slab. A slice reads only the intersected layer of each touched block. Derived fields can require several source arrays per slab. Bin storage scales with the number of bins, and plane storage with the number of transverse pixels. Weighted projections use compensated numerator and denominator maps. No product requires a full three-dimensional field. Text histories and input decks in the separate end-to-end workflow are read independently.
