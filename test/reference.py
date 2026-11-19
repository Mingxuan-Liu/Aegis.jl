"""Independent bounded h5py oracle for the specified private integration suite.

Runs on one complete, unrefined Cartesian snapshot at a time. Reads NumPy
arrays in the writer's (variable, block, z, y, x) order. Emits TOML for Julia.
Optionally validates its samples against a supplied Athena++ reference reader.
"""
import argparse
import importlib.util
import json
import math

import h5py
import numpy as np


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def check_upstream(reader_path, snapshot, field_names, bounds, sample):
    spec = importlib.util.spec_from_file_location("athena_reference", reader_path)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    # Pick an interior point strictly between the faces. For these unrefined
    # outputs the merged reference returns exactly one native cell per field.
    kwargs = {f"x{d}_{edge}": (left + right) / 2
              for d, (left, right) in enumerate(bounds, 1)
              for edge in ("min", "max")}
    merged = reader.athdf(snapshot, quantities=field_names, **kwargs)
    for name, expected in zip(field_names, sample):
        assert merged[name].size == 1, name
        assert float(merged[name].item()) == expected, name


def summarize(path, quantity, reader):
    with h5py.File(path) as h:
        names = [text(n) for n in h.attrs["VariableNames"]]
        datasets = [text(n) for n in h.attrs["DatasetNames"]]
        sizes = h.attrs["NumVariables"]
        assert len(set(names)) == len(names) == sum(sizes)
        fields = [(d, i) for d, n in zip(datasets, sizes) for i in range(n)]
        layout = dict(zip(names, fields))
        nb = int(h.attrs["NumMeshBlocks"])
        size = list(map(int, h.attrs["MeshBlockSize"]))
        sample_block = nb // 3
        indices = [min(n - 1, i) for n, i in zip(size, (5, 11, 19))]
        i, j, k = indices
        samples = [float(h[d][v, sample_block, k, j, i]) for d, v in fields]
        bounds = [list(map(float, h[f"x{d}f"][sample_block, ix:ix + 2]))
                  for d, ix in enumerate(indices, 1)]
        centers = [float(h[f"x{d}v"][sample_block, ix])
                   for d, ix in enumerate(indices, 1)]
        assert text(h.attrs["Coordinates"]) == "cartesian"
        assert np.all(h["Levels"][:] == 0)
        assert np.all(h.attrs["RootGridSize"] % h.attrs["MeshBlockSize"] == 0)
        shape = h.attrs["RootGridSize"] // h.attrs["MeshBlockSize"]
        expected = set(np.ndindex(tuple(shape)))
        assert {tuple(x) for x in h["LogicalLocations"][:]} == expected
        assert nb == len(expected)
        time = float(h.attrs["Time"])
        cycle = int(h.attrs["NumCycles"])
        density_d, density_v = layout["rho"]
        quantity_d, quantity_v = layout[quantity]
        # A small asymmetric box chosen by the native faces of an interior block.
        stops = [min(ix + span, n) for ix, span, n in zip(indices, (2, 3, 4), size)]
        roi = []
        widths = []
        for d, (start, stop) in enumerate(zip(indices, stops), 1):
            faces = h[f"x{d}f"][sample_block, start:stop + 1].astype(np.float64)
            roi.append([float(faces[0]), float(faces[-1])])
            widths.append(np.diff(faces))
        dx, dy, dz = widths
        roi_dv = dz[:, None, None] * dy[None, :, None] * dx[None, None, :]
        i1, j1, k1 = stops
        roi_q = h[quantity_d][quantity_v, sample_block, k:k1, j:j1, i:i1]
        mass, volume, mass_q = [], [], []
        qmin, qmax = math.inf, -math.inf
        for b in range(nb):
            dx, dy, dz = [np.diff(h[f"x{d}f"][b].astype(np.float64))
                          for d in (1, 2, 3)]
            dv = dz[:, None, None] * dy[None, :, None] * dx[None, None, :]
            rho = h[density_d][density_v, b].astype(np.float64)
            q = h[quantity_d][quantity_v, b].astype(np.float64)
            dm = rho * dv
            mass.append(float(dm.sum()))
            mass_q.append(float((dm * q).sum()))
            volume.append(float(dv.sum()))
            qmin = min(qmin, float(q.min()))
            qmax = max(qmax, float(q.max()))

    if reader:
        check_upstream(reader, path, names, bounds, samples)
    # Strings from the supplied files are ASCII and JSON-quoted strings are
    # valid TOML strings. Use finite values only.
    table = dict(names=names, datasets=[d for d, _ in fields], time=time,
                 cycle=cycle, block=sample_block + 1, blocksize=size,
                 indices=[ix + 1 for ix in indices], samples=samples,
                 faces=bounds, centers=centers, volume=math.fsum(volume),
                 mass=math.fsum(mass), mass_mean=math.fsum(mass_q)/math.fsum(mass),
                 minimum=qmin, maximum=qmax, roi=roi, roi_count=int(roi_q.size),
                 roi_volume=float(roi_dv.sum()), roi_integral=float((roi_dv * roi_q).sum()))
    for key, value in table.items():
        print(f"{key} = {json.dumps(value, allow_nan=False)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot")
    parser.add_argument("quantity")
    parser.add_argument("--athena-reader", help="path to compatible athena_read.py")
    args = parser.parse_args()
    summarize(args.snapshot, args.quantity, args.athena_reader)
