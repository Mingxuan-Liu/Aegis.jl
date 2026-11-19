# Newtonian fluid quantities

Declare the gas model at the call site. `NewtonianMHD()` selects stationary Cartesian axes, Newtonian velocity components in one simulation frame, and Athena++ magnetic normalization with permeability one. This declaration is appropriate for the fluid in a Newtonian MHD-PIC calculation even if its particles are relativistic. Input files and field names alone cannot certify the compiled fluid model. Use only primitive density, gas pressure, physical velocity components and cell-centered magnetic components in consistent code units. The default selectors match the inspected writers. Supply other verified selectors explicitly for a custom layout.

| Quantity symbol | Local value | Required inputs |
| --- | --- | --- |
| `speed`, `magnetic_magnitude` | Euclidean norm | Three velocity or three magnetic components |
| `kinetic_energy_density` | $\rho \lVert v\rVert^2/2$ | Density, velocity |
| `magnetic_energy_density`, `magnetic_pressure` | $\lVert B\rVert^2/2$ | Magnetic field |
| `gas_pressure` | $P$, or $\rho c_i^2$ for isothermal gas | Stored pressure, or density and isothermal closure |
| `plasma_beta` | $P/(\lVert B\rVert^2/2)$ | Gas pressure as above, magnetic field |
| `alfven_speed`, `alfven_mach` | $v_A=\lVert B\rVert/\sqrt\rho$, $\lVert v\rVert/v_A$ | Density, magnetic field, velocity for the Mach number |
| `thermal_energy_density` | $P/(\gamma-1)$ | Ideal adiabatic EOS, stored pressure |
| `sound_speed`, `sonic_mach` | $c_s=\sqrt{\gamma P/\rho}$ (adiabatic) or $c_i$ (isothermal), $\lVert v\rVert/c_s$ | Explicit EOS and its required inputs, velocity for the Mach number |
| `temperature` | $(P/\rho)\,v_0^2 m/k_B$ in kelvin | Explicit EOS, mean thermal particle mass $m$ in grams and fluid velocity unit $v_0$ in cm/s; $P/\rho=c_i^2$ for an isothermal EOS |

Pressure and energy density share the code unit of $\rho v^2$. Velocities use the simulation velocity unit. Beta and Mach numbers are dimensionless. A code-unit density or thermal-particle mass conversion cannot be inferred from the particle module's artificial light speed. The model's temperature scale uses the ideal-gas law and $k_B=1.380649\times10^{-16}$ erg/K. An isothermal gas with fixed composition has constant temperature, while stored fluid energy and an adiabatic internal-energy reservoir are outside that closure.

Each native slab is computed in Float64 from only the declared dependencies. A block with $P$ selected cells and source element sizes $s_1,\ldots,s_m$ requires $P(8+\sum_j s_j)$ bytes of input/output array payload, plus HDF5 and Julia bookkeeping. Time is proportional to the selected cell count for each fixed diagnostic. With threaded summaries, the file-owning task evaluates the derived slab before handing the array to an arithmetic worker. The result is local to a cell. `summarize(mesh, q).integral` integrates that local value over the selected native cells. For a nonlinear quantity, the order of operations matters. The methods require positive density where used, nonnegative gas pressure and finite inputs and results. A zero magnetic field yields zero magnetic pressure and Alfvén speed. A zero denominator in beta or either Mach number raises `DomainError` at evaluation, including zero over zero. A nonfinite intermediate or a nonzero ratio or energy that underflows to zero is rejected, even when rescaled algebra could represent the final result. Subnormal values have reduced relative precision.

The magnetic convention and Newtonian energy decomposition follow Stone et al. (2020), §3.1, Eqs. (6)–(8), and Sun & Bai (2023), §§2.1–2.2, Eqs. (2)–(6). The local compatible writer maps the primitive pressure, velocities and cell-centered magnetic array in `OutputType::LoadOutputData`. The EOS classes implement the adiabatic and isothermal sound speeds in `EquationOfState::SoundSpeed`. [Stone et al.](https://doi.org/10.3847/1538-4365/ab929b) · [Sun & Bai](https://arxiv.org/abs/2304.10568) · [Output writer](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/outputs.cpp) · [Adiabatic EOS](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/eos/adiabatic_mhd.cpp) · [Isothermal EOS](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/eos/isothermal_mhd.cpp)


Run `julia --project=. examples/mhd_workflow.jl --demo /tmp/aegis-mhd-example` with a fresh output directory for a public example with known kinetic and magnetic energies. To inspect a thermal diagnostic in a run with a verified adiabatic EOS, the same block interface accepts its EOS parameter explicitly.

```julia
input = read_athinput(input_path)
gas = NewtonianMHD(eos=IdealGas(parameter(input, "hydro", "gamma", Float64)))
open_snapshot(snapshot_path) do s
    cs = quantity(s, gas, :sound_speed)
    read_field(s, cs, 1; ranges=(1:1, 1:1, 1:1))
end
```
