# Athena++ SRMHD writer fixture

The two `.athdf` files are unmodified first outputs of the same special-relativistic Athena++ run. They sample the eight active cells in two Cartesian MeshBlocks at time and cycle zero. The solver used double precision for evolution and HDF5 output with no particles, general relativity, or MPI. The standalone check uses the Aegis public reader and matches primitive and conserved output on the stored coordinates. Run this from the Aegis checkout after instantiating the Julia project.

```sh
julia --project=. test/sr_writer/check.jl
```

## What should the writer store?

`athinput.oblique` gives identical left and right states to the `gr_shock_tube` generator. Its input velocity is the laboratory three-velocity $(0.36,0.48,0)$ in units of $c$. Its speed is $0.6$, so $W=1.25$ and the writer's `vel1`, `vel2`, `vel3` should be $(0.45,0.60,0)$. The rest density and gas pressure are 4 and 1. The independently stored conserved density `dens` should be $D=\rho W=5$. The lab field $(3,4,5)$ has magnitudes 5 both along and perpendicular to the velocity. In the comoving frame these are 5 and $5/W=4$, giving $p'_B=(25+16)/2=20.5$ in code pressure units ($c=\mu=1$). The run includes the producer's initialized primitive recovery, so even a nominally zero velocity component may be of order machine roundoff. This fixture tests the producer's conventions at a single initial state. [Problem generator](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/pgen/gr_shock_tube.cpp#L158-L229) · [Conversion to conserved variables](https://github.com/Mingxuan-Liu/MHD-PIC-in-Athena-/blob/0e3e7a272e1f149d8416020b5ea134bd772242e1/src/eos/adiabatic_mhd_sr.cpp#L296-L325) · [Diagnostic definitions](../../docs/relativity.md)

## Regenerating the snapshot pair

Point `ATHENA_REPO` to a clone of `Mingxuan-Liu/MHD-PIC-in-Athena-` containing commit `0e3e7a272e1f149d8416020b5ea134bd772242e1`. Set `HDF5_PREFIX` to an HDF5 installation with headers in `include` and libraries in `lib`. Run these commands from the Aegis checkout. Export only the source and configuration templates to a fresh temporary directory so that the build uses the selected SR EOS and writer. The supplied build used Apple Clang 21 via the `mpicxx` wrapper for Homebrew HDF5 2.2.0. For a serial HDF5 installation with Clang, set `ATHENA_CXX=clang++`.

```sh
input="$PWD/test/sr_writer/athinput.oblique"
work=$(mktemp -d)
git -C "$ATHENA_REPO" archive 0e3e7a272e1f149d8416020b5ea134bd772242e1 \
    configure.py Makefile.in src | tar -xf - -C "$work"
cd "$work"
python3 configure.py --prob gr_shock_tube --coord cartesian --eos adiabatic \
    -s -b -hdf5 -h5double --cxx clang++ --ccmd "${ATHENA_CXX:-mpicxx}" \
    --hdf5_path "$HDF5_PREFIX"
make -j4
mkdir output
bin/athena -i "$input" -d output
```

Use `output/oblique.prim.00000.athdf` and `output/oblique.cons.00000.athdf`. The `00001` files are the writer's final outputs at the same epoch. Across HDF5 versions, compare variable values and metadata because the file encoding can differ.
