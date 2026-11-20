# Special-relativistic fluid diagnostics

## How fast is the fluid moving relative to the lab?

Run the examples in order from the repository root with `julia --project=.`. This uniform snapshot comes from the SRMHD solver. Its input specifies the lab three-velocity $(0.36,0.48,0)$, so the speed is $0.6c$ and the Lorentz factor should be $1.25$. [Fixture and build instructions](../test/sr_writer/README.md)

```julia
using Aegis

prim = "test/sr_writer/oblique.prim.00000.athdf"
gas = SpecialRelativisticMHD()
motion = open_snapshot(prim) do s
    firstcell(name) = only(read_field(s, quantity(s, gas, name), 1; ranges=(1:1, 1:1, 1:1)))
    (; W=firstcell(:lorentz_factor),
       v=(firstcell(:lab_velocity_1), firstcell(:lab_velocity_2), firstcell(:lab_velocity_3)))
end
@assert motion.W ≈ 1.25
@assert all(isapprox.(motion.v, (0.36, 0.48, 0.0); atol=2e-14))
```

Here $u^\mu=(W,\mathbf u)=W(1,\mathbf v)$ is the four-velocity. In Cartesian coordinates with metric signature $(-,+,+,+)$, $u_\mu u^\mu=-1$. Thus the stored spatial components determine the Lorentz factor and the lab three-velocity. On each output cell $a$ the code evaluates

$$
W_a=\sqrt{1+\sum_{i=1}^3 u_{i,a}^2},
\qquad v_{i,a}=u_{i,a}/W_a.
$$

`:lorentz_factor` is dimensionless. `:lab_velocity_1`, `:lab_velocity_2` and `:lab_velocity_3` return $v_i/c$ along the three simulation axes. A spatial four-velocity component may exceed one. At very high boosts a subluminal three-velocity component can round to one in floating-point arithmetic. Use the stored $\mathbf u$ to calculate $W$ directly. [Stone et al. (2020), §4.1](https://arxiv.org/html/2005.06651v1#S4.SS1) · [Mignone & Bodo (2006), Eqs. (4)–(5)](https://doi.org/10.1111/j.1365-2966.2006.10162.x)

## What magnetic pressure does a fluid element measure?

```julia
pB = open_snapshot(prim) do s
    pressure = quantity(s, gas, :comoving_magnetic_pressure)
    only(read_field(s, pressure, 1; ranges=(1:1, 1:1, 1:1)))
end
@assert pB ≈ 20.5
```

The solver input specifies $\mathbf B=(3,4,5)$ in lab code units with magnetic permeability one. Its components parallel and perpendicular to $\mathbf v$ each have magnitude 5. In the ideal-fluid rest frame these magnitudes become $5$ and $5/W=4$, giving a magnetic pressure of $(5^2+4^2)/2=20.5$ in code energy per unit volume.

The ideal SRMHD magnetic four-vector obeys $b^0=W(\mathbf v\cdot\mathbf B)$, $\mathbf b=\mathbf B/W+b^0\mathbf v$, and $u_\mu b^\mu=0$. The Lorentz scalar $b_\mu b^\mu$ equals the squared rest-frame magnetic field. Using the lab-frame cell-centered magnetic field and spatial four-velocity, the code evaluates the positive sum

$$
p'_{B,a}=\frac{1}{2}b_\mu b^\mu
=\frac{1}{2}\left[\sum_{i=1}^3\left(\frac{B_{i,a}}{W_a}\right)^2+
\left(\frac{\sum_{i=1}^3u_{i,a}B_{i,a}}{W_a}\right)^2\right].
$$

The dot product uses a common power-of-two scale on the velocity components and compensated arithmetic to preserve cancellation for a nearly transverse field. A zero field gives zero pressure. The diagnostic reconstructs the pressure from the saved cell values. Subcell field–velocity correlations and any pressure computed in situ with different sample placement have a different averaging prescription. [Mignone & Bodo (2006), §2, Eqs. (4)–(6)](https://iris.unito.it/retrieve/e27ce42f-1428-2581-e053-d805fe0acbaa/2006.MignoneBodo.MNRAS.pdf)

## Which writer produced these numbers?

`SpecialRelativisticMHD()` is the caller's declaration that the gas primitives follow the flat-spacetime, ideal SR convention in a single Cartesian inertial frame with $c=1$. Check the fluid build, write routine and magnetic normalization for your own file. The compatible writer can transform velocities for orbital advection. Here the verified branch saves ordinary, untransformed primitive slots. GR output, resistive relativistic fluids with an independent electric field, and other velocity conventions require different assumptions. A particle module's artificial light speed characterizes the particles. Gas in a Newtonian MHD-PIC output uses `NewtonianMHD()`. Raw `read_field` takes a stored selector without either model declaration. [Stone et al. (2020), §4.1 and §4.2.1](https://arxiv.org/html/2005.06651v1#S4.SS1) · [Sun & Bai (2023), §§2.1–2.2](https://arxiv.org/abs/2304.10568)

| Native group and field | Physical quantity in this producer | Code units |
| --- | --- | --- |
| `prim/rho`, `prim/press` | Proper density $\rho$ and scalar rest-frame gas pressure $p_g$ | Mass per proper volume and energy per volume |
| `prim/vel1` through `prim/vel3` | Spatial four-velocity $u_i=Wv_i$ in the lab | Dimensionless, $c=1$ |
| `B/Bcc1` through `B/Bcc3` | Cell-centered lab field $B_i$ | Code magnetic field, $B^2$ has energy-density units |
| `cons/dens` | Conserved lab density $D=W\rho$ from the paired conserved output | Mass per lab volume |

The pinned writer copies these scalars and vectors from the hydro and field arrays. Its adiabatic SR EOS computes $W$ from those primitive velocity slots and writes $D=W\rho$. In this example $\rho=4$ and the independently written conserved field gives $D=5$ in every cell. [Primitive and conserved density](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/outputs.cpp#L347-L365) · [Gas pressure](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/outputs.cpp#L387-L395) · [Velocity](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/outputs.cpp#L473-L487) · [Magnetic slots](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/outputs/outputs.cpp#L618-L630) · [SR EOS](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/eos/adiabatic_mhd_sr.cpp#L296-L325) · [Mignone & Bodo (2006), Eq. (10)](https://doi.org/10.1111/j.1365-2966.2006.10162.x)

## Independent checks

The two files in the example are matched primitive and conserved outputs from the pinned SR solver. The public analytic tests separately write manufactured Float32 and Float64 states, transform the parallel and transverse lab fields to the rest frame, and check $u_\mu u^\mu=-1$ and $u_\mu b^\mu=0$. One manufactured state uses the published lab three-velocity $(0.999,0,0)$ and field $(10,7,7)$ from Mignone & Bodo, Table 1, Test 4L. The fixture converts that three-velocity to spatial four-velocity before writing it. Its independent initial-state prediction is $p'_B=[100+98(1-0.999^2)]/2=50.097951$. [Analytic and reader checks](../test/relativistic.jl) · [Mignone & Bodo, Table 1](https://iris.unito.it/retrieve/e27ce42f-1428-2581-e053-d805fe0acbaa/2006.MignoneBodo.MNRAS.pdf)

## Custom outputs and numerical scope

Use `four_velocity=(...)` and, for pressure, `magnetic=(...)` to map three distinct cell-centered fields per vector in x1, x2, x3 order. If your writer saved the equivalent result, `quantity(s, gas, :comoving_magnetic_pressure; stored="my_pb_rest")` declares that field has the required frame, units and centering. The two routes set `q.source` to `:reconstructed` and `:stored`, respectively. Raw identities appear in `q.inputs`. To preserve custom labels, pass selectors returned by `interpret_field`. Derived product metadata carry the labels, model, frame, and selected raw fields. A declaration applies to already compatible numbers and performs no unit conversion.

All operations act independently on corresponding cell indices in the requested native slab. They require no neighboring cells, weighting, periodicity, boundary conditions or full-grid allocation. For $N$ selected cells the working arrays occupy $N(8+\sum_i s_i)$ bytes plus I/O bookkeeping, where $s_i$ is each required input's on-disk element size in bytes. Reconstruction reads three components for $W$ or a velocity and six for pressure. Reading a stored equivalent requires one field. Reconstruction and range checks use Float64. Nonfinite inputs, overflows and nonzero quantities that underflow to zero raise an error. Subnormal Float64 results have reduced precision. Any subsequent spatial averaging must specify its measure.
