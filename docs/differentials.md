# Reconstructed Cartesian derivatives

`DifferentialGrid` builds a lookup from physical MeshBlock location to file block. It takes a verified `CartesianMesh` and three explicit boundary choices. The input and reconstructed output are cell samples in the stored Cartesian frame. A derivative has units of the field divided by one coordinate length. Vector components must be expressed in the same Cartesian basis and compatible units.

- `:periodic` wraps neighboring samples at both ends of the axis. Use it only for a physically periodic field in this frame. Shearing-periodic boundaries require a different remapping.
- `:onesided` reconstructs the derivative at each end from the first or last three interior centers. It can be used at a nonperiodic wall or open boundary when the desired quantity is the interior derivative. This rule imposes no wall parity or ghost state.
- `:invariant` declares zero physical derivative along a singleton direction. A singleton measurement alone cannot establish this symmetry. Every other axis needs at least three distinct centers. The declared invariant contribution is identically zero and reads no samples. The active-axis stencil rejects each nonfinite sample it uses.

For each requested center, fit the quadratic interpolating the field values at three neighboring centers and differentiate it there. If consecutive spacings are $h_1,h_2$, the two polynomial secant slopes are $s_1=(f_b-f_a)/h_1$ and $s_2=(f_c-f_b)/h_2$. The derivative at $x$ is $s_1+(s_2-s_1)[2(x-x_b)+h_1]/(h_1+h_2)$. The stencil is exact for quadratic point samples, and its truncation error for smooth point samples is second order in the maximum local spacing. On a finite-volume mesh the archived value can itself be a cell average. Abrupt width changes or discontinuities can change the convergence to the continuum derivative. The implementation first converts saved samples to Float64 and uses the actual center spacings, including the wrapping interval. [Fornberg (1988)](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

For a vector $u$, `vector_gradient(grid,u)[i][j]` gives $\partial_j u_i$. `divergence(grid,u)` is the tensor trace. `curl(grid,u)[i]` is $\partial_j u_k-\partial_k u_j$ for cyclic $(i,j,k)$. In a Cartesian right-handed frame vorticity is the curl of velocity. A gradient tensor of archived `Bcc` components reconstructs magnetic gradients. In Athena++ constrained transport, native magnetic divergence uses outward-oriented magnetic fluxes on the cell faces divided by cell volume. Computing the trace from `Bcc` samples instead yields a derivative of the archived cell field. It can be nonzero for a solver state satisfying discrete flux conservation. [Stone et al. (2020), §3.2](https://doi.org/10.3847/1538-4365/ab929b) · [Solver conversion to cell-centered fields](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/field/field.cpp#L109) · [Input contract](input-contract.md)

```julia
open_snapshot(path) do s
    mesh = CartesianMesh(s)
    grid = DifferentialGrid(mesh; boundary=(:periodic, :periodic, :periodic))
    omega = curl(grid, ("vel1", "vel2", "vel3"))
    gradB = vector_gradient(grid, ("Bcc1", "Bcc2", "Bcc3"))
    omega3 = read_field(s, omega[3], 1)
    stats = summarize(mesh, gradB[1][2]; region=Region(x1=(0,1)))
end
```

Each directional calculation reads the requested block slab and the necessary neighbor planes, including neighbors outside a selection. For $P$ selected cells and a transverse plane of $A$ cells it uses $O(P+A)$ memory for each term. An edge stencil can add up to two planes. Logical block lookup costs $O(N_{\rm block})$ and saved axis coordinates cost $O(N_1+N_2+N_3)$. Reuse the grid across block reads. Applying successive derivatives to a completed derivative estimates a new operator and widens its effective stencil. The numerical contract here covers first derivatives of the supplied cell samples.

## Convergence and simulation check

`test/differentials.jl` samples $v_y=2.25\sin(2\pi x)$ on the unit periodic interval and compares its reconstructed curl with $\omega_z=4.5\pi\cos(2\pi x)$. It repeats the calculation on one block and four shuffled blocks. The table gives the maximum absolute error at ordinary interior centers, centers touching a block seam, and the two centers touching a domain boundary.

| Axis cells | Interior | Block seams | Periodic domain boundary |
| ---: | ---: | ---: | ---: |
| 24 | 0.14869 | 0.15956 | 0.15956 |
| 48 | 0.03956 | 0.04025 | 0.04025 |
| 96 | 0.01004 | 0.01009 | 0.01009 |

Separate irregular three-dimensional tests exercise nonperiodic boundaries, quadratic samples and blocks as thin as one cell. The optional real check reads a cell at a block seam in the last supplied 128³ filament snapshot and independently differences its upstream and downstream velocity samples. At code time 10 it recovers $\omega_3=0.8564044982$ in inverse code time. Run it with `julia --project=. test/differentials_integration.jl`.
