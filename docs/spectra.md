# Three-dimensional shell spectra

`shell_spectrum` measures the isotropic Fourier spectrum of a complete, uniformly sampled periodic box. Magnetic energy in Athena++ Newtonian code units is $|\mathbf B|^2/2$. Pass the three cell-centered components from the same snapshot. Confirm the physical boundary conditions from the run configuration. The regular sampling test uses the faces and midpoints saved in the snapshot. [Stone et al. (2020), §3.1](https://doi.org/10.3847/1538-4365/ab929b)

```julia
open_snapshot(path) do s
    mesh = CartesianMesh(s)
    periodic = (true, true, true)  # verified against both boundaries on each axis
    B = ("Bcc1", "Bcc2", "Bcc3")
    result = shell_spectrum(mesh, B; periodic)
    sum(result.shell_energy) + result.dc_energy + result.excluded_energy
end
```

For $N=N_1N_2N_3$ samples, normalize each three-dimensional forward coefficient by $N$. Each shell stores

$$
S_j = \tfrac12 \sum_{a} \sum_{\mathbf{k}\in j}|\widehat f_a(\mathbf{k})|^2,
\qquad M_j=S_j/\Delta k_j,\qquad k_d=2\pi m_d/L_d .
$$

The caller must establish a common unit and physical frame for the components $a$, including unlabeled custom fields. Use a single-element tuple for a scalar. With `density="rho"` and three Newtonian velocities the transformed fields are $\sqrt{\rho}\,v_a$. Supplying this selector declares that the field is gas mass density in compatible units. Their integrated spectral energy measures $\langle\rho |\mathbf v|^2/2\rangle$. Subtracting the mean in this case subtracts each transformed field's spatial mean. Multiplying a velocity spectrum by the mean density is a different diagnostic. Labels can record the user's units, and distinct supplied component unit labels raise an error.

| Returned value | Convention |
| --- | --- |
| `edges`, `k` | Angular wavenumber in radians per code length. Default bins are linear to the Fourier-grid corner. Callers may supply nonnegative edges and choose arithmetic or geometric centers. |
| `shell_energy`, `density` | $S_j$ and $M_j$. A shell includes its left edge. The last shell also includes its right edge. An empty shell has zero energy, density and modes. |
| `mode_counts`, `excluded_modes` | Integer counts of full-grid modes inside shells and outside the supplied edges. DC is a separate single mode. Counts are independent of the number of components. |
| `dc_energy`, `original_mean_energy` | Measured DC of the transformed field and energy in the original component means. Setting `remove_mean=true` subtracts each component mean in real space before its FFT. |
| `real_energy`, `excluded_energy` | Independent real-space half squared norm per cell of the transformed fields and the summed energy of non-DC modes outside the edges. |
| `provenance` | Snapshot path/time/cycle, field and optional density identities, grid dimensions and lengths, periodicity declaration, FFT thread count and workspace estimate. |

For components with unit $U$, shell totals and real-space energy have units $U^2$. A density per angular wavenumber has units $U^2$ times code length. For density-weighted components both expressions gain a factor of density units. The magnetic reference figure excludes DC from the displayed shells while retaining the original field for the transform.

The array is assembled by logical block location with physical $x_1$ varying fastest. The global real transform compresses the first Julia dimension. For the full shell power, an interior positive $k_1$ plane carries weight two. The $k_1=0$ plane and the even-length $k_1$ Nyquist plane carry weight one, because those planes retain all $k_2,k_3$. For odd $N_1$ the largest positive $k_1$ also carries weight two. Each mode's normalized squared modulus is accumulated with a compensated shell sum. With all modes accounted for, $\sum_j M_j\Delta k_j + E_{\mathrm{DC}} + E_{\mathrm{excluded}} = E_{\mathrm{real}}$ to rounding precision. [FFTW reference, §4.8.6](https://www.fftw.org/fftw3_doc/Multi_002ddimensional-Transforms.html) · [AbstractFFTs real transforms](https://juliamath.github.io/AbstractFFTs.jl/stable/api/)

For an inscribed spherical cutoff, supply a final edge at the desired angular radius. The default retains the high-wavenumber corners where shell coverage is directionally incomplete. Switching a plotted axis to $q=kL/(2\pi)$ leaves the physical ordinate unchanged if the axis alone changes. To report a density **per unit $q$**, divide shell energy by $\Delta q$, equivalent to $M_q=(2\pi/L)M_k$. State which wavenumber and density enter any compensated ordinate.

## Limits and cost

Prepare a complete unrefined 3D Cartesian mesh. Every axis must have at least two cells. Along each axis, faces and centers must coincide with the uniform lattice within four coordinate-roundoff units on the box scale, capped at 0.002 cell widths. Nonuniform spacing, AMR or a nonperiodic physical field requires another estimator. Archive precision, stale boundary declarations and a mismatch between stored variables and their intended physical model can affect interpretation. Expect rounding for finite Float64 arithmetic, and choose scaled input units when the squared norm exceeds its dynamic range.

The transform is global. One Float64 grid, one ComplexF64 half-grid and one FFTW ESTIMATE plan are reused for each component. The algorithm visits the grid for each transform and performs a shell lookup per retained Fourier coefficient. For $C$ components and $J$ bins, work scales as $O(CN(\log N + \log J))$. Use `spectrum_workspace(s; bins=J)` before requesting a large case. The preflight estimate reserves a further real-grid buffer for FFT workspace and two decoded blocks for optional density. The default budget is 1 GiB. Platform-specific FFTW work buffers, runtime memory, storage caches and derived-field scratch space also affect process memory.

| Grid | Estimated workspace for $32^3$ blocks and 200 shells |
| --- | ---: |
| $128^3$ | 0.048 GiB |
| $256^3$ | 0.376 GiB |
| $512^3$ | 3.00 GiB |
| $1024^3$ | 24.02 GiB |
| $2048^3$ | 192.06 GiB |

The tests compare against a full direct complex DFT on odd/even, oblique, unequal-length and non-cubic lattices, with known amplitudes for interior, DC and Nyquist modes. The [magnetic spectrum comparison](../validation/spectra/) uses the supplied filament-heating outputs.
