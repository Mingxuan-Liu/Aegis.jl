# Three-dimensional shell spectra

How much magnetic or Newtonian kinetic energy lies at a given spatial scale? Run this example from the repository root after installing the dependencies. It creates a periodic $8^3$ cube, computes both spectra and prints their energy closures.

```sh
julia --project=. examples/periodic_spectra.jl
```

The same calls apply to a simulation snapshot once you verify its Cartesian geometry, periodic boundaries and primitive fields. To try them on the manufactured wave in a Julia session started with `--project=.`, run the following from the repository root.

```julia
using Aegis
include("examples/periodic_spectra.jl")
mktempdir() do dir
    path = PeriodicSpectraExample.write_wave(joinpath(dir, "wave.athdf"))
    open_snapshot(path) do s
        mesh = CartesianMesh(s)
        periodic = (true, true, true)
        edges = 1.0:2.0:15.0  # radians per code length in this example
        magnetic = shell_spectrum(mesh, ("Bcc1", "Bcc2", "Bcc3"); periodic, edges)
        kinetic = shell_spectrum(mesh, ("vel1", "vel2", "vel3"); density="rho", periodic, edges)
        q = magnetic.k .* (π / (2π))  # this cube has side L = π
        compensated = magnetic.k.^(5/3) .* magnetic.density
        (; first_q=first(q), first_ordinate=first(compensated),
           magnetic_energy=magnetic.real_energy, kinetic_energy=kinetic.real_energy)
    end
end
```

The energies in the simulation frame are $\langle |\mathbf B|^2/2\rangle_V$ and $\langle\rho |\mathbf v|^2/2\rangle_V$, where $\rho$ is the gas mass density, $\mathbf v$ is the gas three-velocity, and the magnetic permeability is one in code units. Both are energy per volume. Angular brackets denote the spatial mean of the pointwise energy. On the uniform lattice this is an equal-cell average. [Stone et al. (2020), §3.1](https://doi.org/10.3847/1538-4365/ab929b)

For $N=N_1N_2N_3$ cell centers indexed by $\mathbf n$ starting at zero, form one Fourier input $f_a$ per component. Use $f_a=B_a$ for magnetic energy and $f_a=\sqrt{\rho}\,v_a$ for kinetic energy. Given box lengths $L_d$, a signed Fourier index $m_d$ corresponds to angular wavenumber $k_d=2\pi m_d/L_d$. The transform uses the first sample as its phase origin. A common translation multiplies the coefficients by a phase and preserves shell energies.

$$
\widehat f_a(\mathbf m)=\frac{1}{N}\sum_{\mathbf n}f_a(\mathbf n)
\exp\!\left[-2\pi i\sum_{d=1}^3\frac{m_dn_d}{N_d}\right],
\qquad
S_j=\frac{1}{2}\sum_a\sum_{\substack{\mathbf m\ne\mathbf 0\\|\mathbf k|\in I_j}}|\widehat f_a(\mathbf m)|^2,
\qquad M_j=\frac{S_j}{\Delta k_j}.
$$

Here $I_j$ is a radial bin with width $\Delta k_j$. Its integral $S_j$ (`shell_energy`) is energy per volume summed over modes. Its density $M_j$ (`density`) is energy per volume *per unit angular wavenumber*. Spherical integration applies to anisotropic fluctuations too, but the shell totals discard wavevector direction. With all DC and out-of-range power accounted for, $\sum_j M_j\Delta k_j+E_{\rm DC}+E_{\rm excluded}=\langle\sum_a f_a^2/2\rangle_V$ to rounding precision. With centering enabled, use the centered $f_a$ on the right.

### Check a wave by hand

In this cube $L=\pi$ and the fields at each cell center $x$ are

$$
\mathbf B=(3,2\sin(2x),0),\qquad \rho=(2+\sin(2x))^2,\qquad \mathbf v=(0,0,2).
$$

A sinusoid of amplitude 2 has coefficients of magnitude 1 at each of the two modes $k_x=\pm2$. Together they contribute $S_1=(1^2+1^2)/2=1$. In the first bin $[1,3)$, $M_1=1/2$. The filament-paper convention uses $q=kL/(2\pi)$ on the horizontal axis and $k^{5/3}M_j$ on the vertical axis, evaluated at the shell center $k$. Here $k=2$ gives $q=1$ and the compensated ordinate is $2^{2/3}\simeq1.58740105$.

| FFT input | DC energy | First shell | Real-space energy |
| --- | ---: | ---: | ---: |
| Magnetic field $\mathbf B$ | 4.5 | 1 | 5.5 |
| Weighted velocity $\sqrt{\rho}\,\mathbf v=(0,0,4+2\sin(2x))$ | 8 | 1 | 9 |
| Weighted velocity with `remove_mean=true` | 0 | 1 | 1 |
| Gas velocity in the rest frame of this uniform flow, then weighted | 0 | 0 | 0 |

The positive-wavenumber shells always exclude the zero mode. Adding a uniform offset to $\mathbf B$, or to any field already prepared for the FFT, therefore leaves these shells unchanged up to roundoff. Setting `remove_mean=true` centers each *prepared component*. For kinetic energy this removes $\langle\sqrt{\rho}\,v_a\rangle_V$. Subtracting a uniform gas velocity $U_a$ *before* weighting adds $-\sqrt{\rho}\,U_a$ to that component. At variable density this term can carry power at nonzero wavenumber. In the example the gas flows uniformly, so its rest-frame kinetic energy is exactly zero even though centering the weighted field leaves the sinusoid. The example uses `cell_map` to subtract that known gas velocity before the FFT. [Executable example](../examples/periodic_spectra.jl)

### Reading the result

The three component selectors must describe co-located fields in a common unit and physical frame. Passing `density="rho"` declares the density in compatible units. For a scalar spectrum use a single-element tuple. Labels can record units and distinct supplied component unit labels raise an error.

| Returned value | Convention |
| --- | --- |
| `edges`, `k` | Angular wavenumber in radians per code length. Default bins are linear to the Fourier-grid corner. Callers may supply nonnegative edges and choose arithmetic or geometric centers. |
| `shell_energy`, `density` | $S_j$ and $M_j$. A shell includes its left edge. The last shell also includes its right edge. An empty shell has zero energy, density and modes. |
| `mode_counts`, `excluded_modes` | Integer counts of full-grid modes inside shells and outside the supplied edges. DC is a separate single mode. Counts are independent of the number of components. |
| `dc_energy`, `original_mean_energy` | Measured DC of the transformed field and energy in the original component means. Setting `remove_mean=true` subtracts each component mean in real space before its FFT. |
| `real_energy`, `excluded_energy` | Independent real-space half squared norm per cell of the transformed fields and the summed energy of non-DC modes outside the edges. |
| `provenance` | Snapshot path/time/cycle, field and optional density identities, grid dimensions and lengths, periodicity declaration, FFT thread count and workspace estimate. |

For components with unit $U$, shell totals and real-space energy have units $U^2$. A density per angular wavenumber has units $U^2$ times code length. For density-weighted components both expressions gain a factor of density units.

The array is assembled by logical block location with physical $x_1$ varying fastest. The global real transform compresses the first Julia dimension. For the full shell power, an interior positive $k_1$ plane carries weight two. The $k_1=0$ plane and the even-length $k_1$ Nyquist plane carry weight one, because those planes retain all $k_2,k_3$. For odd $N_1$ the largest positive $k_1$ also carries weight two. Each mode's normalized squared modulus is accumulated with a compensated shell sum. [FFTW reference, §4.8.6](https://www.fftw.org/fftw3_doc/Multi_002ddimensional-Transforms.html) · [AbstractFFTs real transforms](https://juliamath.github.io/AbstractFFTs.jl/stable/api/)

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

The tests compare against a full direct complex DFT on odd/even, oblique, unequal-length and non-cubic lattices, with known amplitudes for interior, DC and Nyquist modes. The [spectral tests](../test/spectra.jl) generate their own fixtures.
