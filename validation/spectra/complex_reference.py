"""Check a measured, cubic final snapshot using SciPy's full complex FFT."""
import argparse
import csv
from pathlib import Path

import h5py
import numpy as np
from scipy.fft import fftn, fftfreq


def check(snapshot, table):
    with table.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream)
                if row["aegis_snapshot"] == snapshot.name]
    if not rows:
        raise ValueError(f"No measured shells in {table} for {snapshot.name}")
    edges = np.array([float(row["k_lower"]) for row in rows] + [float(rows[-1]["k_upper"])])
    expected = np.array([float(row["aegis_shell_energy"]) for row in rows])

    with h5py.File(snapshot) as h:
        n, ny, nz = map(int, h.attrs["RootGridSize"])
        if n != ny or n != nz or n > 256 or np.any(h["Levels"][:] != 0):
            raise ValueError("This full-complex check is sized for level-zero cubic grids up to 256³")
        assert h.attrs["Coordinates"].decode() == "cartesian"
        assert float(h.attrs["Time"]) == float(rows[0]["aegis_time"])
        lower = np.array([h.attrs[f"RootGridX{d}"][0] for d in (1, 2, 3)])
        lengths = np.array([h.attrs[f"RootGridX{d}"][1] for d in (1, 2, 3)]) - lower
        assert np.all(lengths == float(rows[0]["L"]))
        dx = float(lengths[0]) / n

        # Place blocks by their physical faces. The independent HDF5 reader
        # and SciPy both work with contiguous arrays in (z, y, x) order.
        starts = [np.rint((h[f"x{d}f"][:, 0] - lower[d-1])/dx).astype(int)
                  for d in (1, 2, 3)]
        bx, by, bz = map(int, h.attrs["MeshBlockSize"])
        slots = [(name.decode(), i) for name, count in
                 zip(h.attrs["DatasetNames"], h.attrs["NumVariables"]) for i in range(count)]
        names = [v.decode() for v in h.attrs["VariableNames"]]
        k = (2*np.pi/lengths[0]) * (fftfreq(n)*n)
        kxy2 = k[:, None]**2 + k[None, :]**2
        energy = np.zeros(len(expected))
        dc, real = 0.0, 0.0
        grid = np.empty((n, n, n), dtype=np.float64)
        for component in ("Bcc1", "Bcc2", "Bcc3"):
            assert names.count(component) == 1
            dataset, slot = slots[names.index(component)]
            grid.fill(np.nan)
            for b, (x, y, z) in enumerate(zip(*starts)):
                grid[z:z+bz, y:y+by, x:x+bx] = h[dataset][slot, b]
            assert np.isfinite(grid).all()
            real += np.einsum("i,i", grid.ravel(), grid.ravel()) / (2*grid.size)
            coeff = fftn(grid, norm="forward", workers=1)
            dc += abs(coeff[0, 0, 0])**2/2
            for iz in range(n):
                magnitude = np.sqrt(kxy2 + k[iz]**2)
                energy += np.histogram(magnitude, bins=edges, weights=abs(coeff[iz])**2/2)[0]
            del coeff
    # Histograms of weighted powers and block sums use different summation
    # orders. Small tail bins have the largest relative rounding error.
    np.testing.assert_allclose(dc + energy.sum(), real, rtol=5e-11, atol=0)
    np.testing.assert_allclose(energy, expected, rtol=5e-11, atol=0)
    print(f"{snapshot.name}: maximum relative shell difference {max(abs(energy/expected-1)):.3g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshots", type=Path, nargs="+")
    parser.add_argument("--table", type=Path, default=Path(__file__).with_name("magnetic_shells.csv"))
    args = parser.parse_args()
    for snapshot in args.snapshots:
        check(snapshot, args.table)
