# Magnetic fluctuation energy along and across a guide field

How much magnetic energy has accumulated in fluctuations transverse to the initial guide field? Evaluate the energy in each native cell, then average it over the desired volume. From the repository root, run the following in Julia started with `julia --project=.`.

```julia
using Aegis
include("examples/mhd_workflow.jl")
mktempdir() do dir
    open_snapshot(MHDWorkflowExample.write_demo(dir)) do s
        model = NewtonianMHD()
        background = (3, 0, 0)
        along = magnetic_fluctuation(s, model, :parallel; reference=background)
        across = magnetic_fluctuation(s, model, :transverse; reference=background)
        mesh = CartesianMesh(s)
        @assert summarize(mesh, along).mean == 4.5
        @assert summarize(mesh, across).mean == 2
    end
end
```

Here every cell contains $\mathbf B=(0,0,2)$, so subtracting $\mathbf B_0=(3,0,0)$ gives $\delta\mathbf B=(-3,0,2)$. Along the $x$ axis the signed perturbation is $-3$, giving energy density $9/2$. Across it, the perturbation has magnitude 2, giving energy density 2. `direction` defaults to the given reference vector. To measure the total field relative to an independent axis, use, for example, `reference=(0,0,0), direction=(1,0,0)`. Neither parameter is estimated from the snapshot.

For a uniform reference $\mathbf B_0$ in code magnetic units and a nonzero Cartesian direction $\mathbf d$, define $\hat{\mathbf n}=\mathbf d/|\mathbf d|$. In the Newtonian simulation frame, with Athena++ magnetic permeability one, the parallel and transverse fluctuation energy densities are

$$
e_{\parallel}=\tfrac12(\delta\mathbf B\cdot\hat{\mathbf n})^2,\qquad
e_{\perp}=\tfrac12|\delta\mathbf B\times\hat{\mathbf n}|^2,\qquad
\delta\mathbf B=\mathbf B-\mathbf B_0.
$$

At each cell $i$ the code substitutes the three stored cell-centered components $B_{i,1},B_{i,2},B_{i,3}$ in these equations. The output field is named `magnetic_fluctuation_parallel` or `magnetic_fluctuation_transverse` in code energy density units ($\rho v^2$). No neighbor, ghost sample, boundary condition or Fourier transform enters. The volume average is $\sum_i e_i\Delta V_i/\sum_i\Delta V_i$. Squaring a field averaged within a cell can give a different result from averaging a squared field over that cell. For a circularly polarized perturbation with amplitude $A$ entirely transverse to the guide, $e_\perp=A^2/2$ at every phase. Following this quantity through time measures fluctuation energy growth or damping. The scalar energy discards polarization handedness. [Bell (2004), §2, Eq. (3)](https://adsabs.harvard.edu/pdf/2004MNRAS.353..550B) · [Stone et al. (2020), §3.1, Eq. (8)](https://arxiv.org/html/2005.06651v1#S3.SS1)

By default the reconstruction reads `("Bcc1","Bcc2","Bcc3")`. A verified custom component layout can supply `magnetic=(bx,by,bz)`. These selectors must resolve to distinct stored cell-centered components in the same units and Cartesian frame. The Athena++ output writer obtains these components from `pfld->bcc`. Raw names, dataset slots and any supplied interpretation labels are preserved. [Athena++ compatible MHD writer](https://github.com/dfielding14/athena-1/blob/d2aa0b0ffc3b6b215721693cf82946cfa5582616/src/outputs/outputs.cpp#L1052-L1100)

If the producer already saved an equivalent scalar, pass, for example, `stored=interpret_field(s, "bell_Eperp"; description="pgen transverse fluctuation energy")` to the transverse call. This explicitly declares the same reference, direction, units, centering and snapshot frame. Verify the pgen's definition, normalization and averaging in its source. Aegis reads just that nonnegative scalar and records the source as `stored`. For reconstructed results it records `reconstructed` and the original component selectors. These declarations accompany the field metadata of exported slices, projections and distributions.

For $P$ selected cells, reconstructing one energy reads three input slabs and writes one Float64 slab in $O(P)$ time. A stored equivalent reads one input slab. The transverse calculation evaluates a cross product directly to retain small perpendicular fluctuations in a strong parallel field. Finite precision of the stored components limits the smallest measurable perturbation. Input components that already rounded to the background cannot recover that signal.

On the supplied 256³ MHD snapshot with 512 native 32³ blocks and Float32 magnetic components, a warmed transverse `read_field(s, across, 1)` allocated 661,600 bytes, including its 655,360-byte input/output array payload. A serial full-domain `summarize` took 0.325 s on an Apple M2 Max with Julia 1.12.0 and a warm filesystem cache, using `reference=(0.25,-0.125,0.0625)` and `direction=(3,4,0)`. The public tests exercise the same block shape and include an allocation bound.
