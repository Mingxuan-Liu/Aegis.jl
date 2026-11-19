# Final-snapshot magnetic spectra

The 256³ snapshot supplied under `test/test_data/3D_Athena` reproduces the reference 256³ curve shell by shell. The supplied 128³ snapshot reaches a higher compensated amplitude than the reference 128³ curve. Its resistivity and viscosity are both 0.005, while the plotted run uses 0.013. Two further raw snapshots in the filament project archive allowed like-for-like comparisons at 128³ and 512³. All four measurements use the last supplied state in the corresponding run directory.

[Overlay of native Aegis shells and the figure curves](magnetic_spectra_overlay.pdf) · [Unsmoothed per-shell measurements and residuals](magnetic_shells.csv)

## Comparison by curve

The ordinate is $C=k^{5/3}M(k)$, with $M(k)$ per unit angular wavenumber. The slope below measures the **uncompensated** spectrum across the same finite band in every run. It equals $\log[C(16)/C(4)]/\log(4)-5/3$, evaluated on the log-log polyline of the reference plotting recipe, where the arguments denote $q=kL/(2\pi)$. Both $q=4$ and $q=16$ are inside the measured range. The spectrum curves through this band, so the secant describes this particular range and epoch.

| $N^3$ and $\eta=\nu$ in figure | Figure output and code time | Aegis / figure ordinate at $q=32$ | Aegis / figure slope, $4\le q\le16$ | Largest relative ordinate difference | Outcome |
| --- | --- | ---: | ---: | ---: | --- |
| 128³, 0.013 | 00010, 10 | 0.09781 / 0.09781 | −1.241 / −1.241 | $1.55\times10^{-15}$ | Matching final snapshot from filament archive |
| 256³, 0.005 | 00010, 10 | 0.25723 / 0.25723 | −1.166 / −1.166 | $1.33\times10^{-15}$ | Matching final snapshot supplied in Aegis |
| 512³, 0.002 | 00010, 10 | 0.34698 / 0.34698 | −1.257 / −1.257 | $4.44\times10^{-16}$ | Matching final snapshot from filament archive |
| 1024³, 0.0008 | 00050, 5 | — / 0.42270 | — / −1.128 | — | Original shell table available, raw snapshot absent |
| 2048³, 0.0003 | 00056, 5.26000023 | — / 0.43773 | — / −1.139 | — | Original shell table available, raw snapshot absent. Its native level is 1. |

The matched differences are reported relative to each original shell value. Matched runs have the same input-deck SHA-256, source basename and byte count, grid, time, box size and transport coefficients as the reference records. The match compares independent transforms with no fitted amplitude. The 128³, 256³ and 512³ figure snapshots have $t/t_{\rm eddy}=4.7390$, the 1024³ case has 2.3695, and the 2048³ case has 2.4927.

For the **additional supplied 128³ run** at $\eta=0.005$ and $t=10$, Aegis gives $C(32)=0.20942$ and a slope of **−1.103** on this band. It shares the broad rise and high-$q$ fall of the reference compensated curve. Its maximum unsmoothed shell ordinate lies at $q=17.77$, compared with $q=10.53$ for the 128³ run at $\eta=0.013$. The latter is reproduced by Aegis from its own raw snapshot. For the reference 128³, 256³ and 512³ cases the unsmoothed maxima shift to $q=10.53,19.51,36.65$. Both resolution and diffusivity change along that sequence.

## Measurement and reference provenance

- `measure.jl` opens the field stream with Aegis, verifies a periodic Cartesian cube, reconstructs native spatial placement, and calls `shell_spectrum` with all three saved cell-centered magnetic components. It uses a global Float64 real FFT and keeps the mean-field energy separate. The nonzero shells contain $N^3-1$ full-grid modes including the corner. Including DC, Parseval closes to within $3\times10^{-16}$ in absolute code energy density in every measured run.
- Original numerical reference: the filament project’s `Reconnection-Analysis/production/cache/highres_magnetic_spectra/{run}.h5`. Each of the five caches has a checked SHA-256 over the numerical arrays. Their recorded estimator digest is `928308e66492097ff7f8f8c8d5d61eb17a6549a9859b242b0eb99f88ccad9d04`. `resolution_magnetic_spectra.jl`, `single_magnetic_spectrum.jl`, `spectrum_density.jl`, `large_magnetic_spectrum.jl` and `smoothie.jl` in the source project establish the pipeline (checkout `1b2db089cff2c0fc7bf69739cdfaa64259ee5e5c`).
- All cases have $L=8$ in code length. Each run uses $\lceil4\log_{10}(k_{\rm end}/k_{\min})\rceil$ logarithmic shells from the fundamental $k_{\min}=2\pi/L$ through $k_{\rm end}=1.00001\sqrt{3}(N/2)k_{\min}$. Geometric shell centers define both $q$ and $k^{5/3}$; the plotted density is $S_j/\Delta k_j$. The source uses natural cubic interpolation of the positive compensated native shell measurements in log-log coordinates, sampled at 128 equally spaced log abscissae for drawing. The CSV retains the native values, bin edges, mode counts, uncompensated and compensated densities, and signed residuals. Differences for the distinct-parameter case compare physical runs and have a separate status in the CSV.
- `Paper/figures/magnetic_spectra_128_256_512_1024_2048.pdf`, SHA-256 `2880f55242d0541e980c9a07bb6d38bc67655f16472e9d4122204788ef4a770e`, is the five-curve figure. Calibrating the vector paths on the four major ticks per logarithmic axis gives $C(32)=0.09780,0.25724,0.34695,0.42269,0.43774$ in resolution order. Agreement to the drawing precision ties the numerical caches to the supplied plot.
- An independent full-complex SciPy FFT on the two supplied final snapshots agrees with the Aegis shell energies to maximum relative differences of $3.03\times10^{-13}$ (128³) and $2.61\times10^{-12}$ (256³).
- With one FFT thread, the 512³ transform of three stored components took 7.06 s in the measured call on the 32 GB M2 Max. The Julia process reached 2.98 GiB peak RSS. The preflight estimate for 512³ is about 3.01 GiB, including a full real-grid allowance for FFTW scratch. The same estimator needs about 24 GiB at 1024³ and 192 GiB at 2048³. These exceed the 6 GiB analysis budget used here. The 2048³ cache also records level-one blocks. Applying the native transform to that raw mesh would require validating its refined cell coverage and placement. Timings depend on filesystem cache and compilation.

## Reproduce

From the Aegis repository root, set `FILAMENT` to the local project directory containing `Data`, `Paper` and `Reconnection-Analysis`. The Python reporting script uses NumPy, SciPy, h5py and Matplotlib. The public package’s FFT needs Julia and its declared dependencies.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
measurements=$(mktemp -d)
julia --project=. validation/spectra/measure.jl "$measurements" \
    test/test_data/3D_Athena/128_eta5e-3 test/test_data/3D_Athena/256_eta5e-3 \
    "$FILAMENT/Data/128_eta1.3e-2" "$FILAMENT/Data/512_eta2e-3"
python3 validation/spectra/compare.py \
    --cache "$FILAMENT/Reconnection-Analysis/production/cache/highres_magnetic_spectra" \
    --measurements "$measurements" --out "$measurements/report"
```

To independently check the measured shells on both supplied 3D runs with a full complex FFT, run the following from the same directory.

```sh
python3 validation/spectra/complex_reference.py test/test_data/3D_Athena/*/*00010.athdf
```

To redraw the saved results from the CSV alone, run the following command.

```sh
python3 validation/spectra/compare.py --table validation/spectra/magnetic_shells.csv --out /tmp/aegis-spectrum-report
```
