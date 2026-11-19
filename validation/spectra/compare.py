"""Compare native Fourier shells to the five figure caches and render their overlay."""
import argparse
import csv
import hashlib
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import CubicSpline


def string(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def assemble(cache_dir, measured_dir):
    """Pair each measured state with its same-resolution figure case."""
    measured = {}
    for path in sorted(measured_dir.glob("*.h5")):
        with h5py.File(path) as h:
            measured.setdefault(h.attrs["n"], []).append(path)
    rows = []
    seen = set()
    for path in sorted(cache_dir.glob("*.h5")):
        with h5py.File(path) as h:
            a = dict(h.attrs)
            q, edges, energy = (h["spectra"][key][:] for key in
                                ("kn", "k_edges", "magnetic_shell_energy"))
        n, L = int(a["nx1"]), float(a["Lx1"])
        if n in seen:
            raise ValueError(f"Duplicate resolution in figure caches: {n}")
        seen.add(n)
        digest = hashlib.sha256(np.concatenate([q, edges, energy]).astype("<f8").tobytes())
        assert digest.hexdigest() == string(a["arrays_sha256"]), path
        np.testing.assert_allclose(energy.sum() + a["dc_energy"], a["real_energy"], rtol=5e-9)
        k = q * (2*np.pi/L)  # Restore angular k from the original cached horizontal axis.
        np.testing.assert_allclose(k, np.sqrt(edges[:-1]*edges[1:]), rtol=1e-13)
        ref_density = energy / np.diff(edges)
        ref_ordinate = energy * (k**(5/3) / np.diff(edges))
        for measurement in measured.get(n, [None]):
            data, attrs = {}, {}
            if measurement is not None:
                with h5py.File(measurement) as h:
                    attrs = dict(h.attrs)
                    data = {key: h[key][:] for key in h}
                np.testing.assert_allclose(data["edges"], edges, rtol=2e-14, atol=0)
                np.testing.assert_allclose(data["k"], k, rtol=2e-14, atol=0)
                assert attrs["L"] == L and attrs["n"] == n
                key_pairs = (("run", "run"), ("source", "source_basename"),
                             ("deck_sha256", "deck_sha256"), ("source_bytes", "source_bytes"),
                             ("t_code", "t_code"), ("eta", "eta_ohm"), ("nu", "nu_iso"))
                matched = all(string(attrs[x]) == string(a[y]) for x, y in key_pairs)
                status = "matched" if matched else "distinct run"
                ordinate = data["k"]**(5/3) * data["density"]
                if matched:
                    np.testing.assert_allclose(data["shell_energy"], energy, rtol=2e-10, atol=0)
                    np.testing.assert_allclose(attrs["dc_energy"], a["dc_energy"], rtol=2e-10)
                    np.testing.assert_allclose(attrs["real_energy"], a["real_energy"], rtol=2e-10)
            else:
                status = "reference only"
            for j in range(len(k)):
                tested = measurement is not None
                rows.append(dict(
                    N=n, status=status, reference_snapshot=string(a["source_basename"]),
                    reference_eta=a["eta_ohm"], reference_time=a["t_code"], L=L,
                    aegis_snapshot=string(attrs["source"]) if tested else "",
                    aegis_eta=attrs.get("eta", ""), aegis_time=attrs.get("t_code", ""),
                    shell=j+1, k_lower=edges[j], k_upper=edges[j+1], k=k[j], q=q[j],
                    reference_shell_energy=energy[j], reference_M=ref_density[j],
                    reference_compensated=ref_ordinate[j],
                    fourier_modes=int(data["mode_counts"][j]) if tested else "",
                    aegis_shell_energy=data["shell_energy"][j] if tested else "",
                    aegis_M=data["density"][j] if tested else "",
                    aegis_compensated=ordinate[j] if tested else "",
                    difference=ordinate[j]-ref_ordinate[j] if tested else "",
                    fractional_difference=ordinate[j]/ref_ordinate[j]-1 if tested else ""))
    if seen != {128, 256, 512, 1024, 2048} or set(measured)-seen:
        raise ValueError("Expected the five figure resolutions and no unpaired measurements")
    if not {128,256}.issubset(measured):
        raise ValueError("The comparison requires both supplied 3D_Athena runs")
    return sorted(rows, key=lambda r:(r["N"], r["status"], r["shell"]))


def curve(q, y):
    """Natural cubic interpolation in (log q, log y) at 128 log-spaced points."""
    grid = np.geomspace(q[0], q[-1], 128)
    spline = CubicSpline(np.log(q), np.log(y), bc_type="natural")
    return grid, np.exp(spline(np.log(grid)))


def interp_line(grid, values, point):
    # Interpolate the actual plotted polyline in log coordinates.
    return np.exp(np.interp(np.log(point), np.log(grid), np.log(values)))


def render(rows, output):
    # Left panel collects the figure's resolution sequence. The right panel
    # separates the additional 128³ run at a different resistivity.
    colors = {128:"#245687", 256:"#b15815", 512:"#247969", 1024:"#817985", 2048:"#a29aa3"}
    fig, (ax, small) = plt.subplots(1, 2, figsize=(9.6, 4.25),
                                    gridspec_kw={"width_ratios":[1.48, 1]}, layout="constrained")
    ax.set_title("Reference resolution sequence", loc="left", fontweight="semibold")
    small.set_title("128³ at two resistivities", loc="left", fontweight="semibold")
    for axis in [ax, small]:
        axis.set(xscale="log", yscale="log", ylim=(0.0017,0.7), xlabel=r"$q=kL/(2\pi)$")
        axis.grid(which="major", color="#c9cbd0", lw=0.55, alpha=0.55)
        axis.spines[["top","right"]].set_visible(False)
    ax.set(xlim=(1,1550), ylabel=r"$k^{5/3}\,S_j/\Delta k_j$  (code units)")
    small.set_xlim(1, 103)
    plotted = set()
    groups = {}
    for row in rows:
        groups.setdefault((int(row["N"]),row["status"],row["aegis_snapshot"]), []).append(row)

    for (n,status,source), group in groups.items():
        group.sort(key=lambda r:int(r["shell"]))
        q = np.array([float(r["q"]) for r in group])
        ref = np.array([float(r["reference_compensated"]) for r in group])
        xr, yr = curve(q,ref)
        if n not in plotted:
            # Unavailable levels are shown as reference curves only.
            ax.plot(xr,yr, color=colors[n], lw=1.55,
                    ls={1024:"--",2048:":"}.get(n,"-"),
                    label=rf"{n}$^3$" + ("  (ref. only)" if status == "reference only" else ""))
            plotted.add(n)
            if n == 128:
                small.plot(xr,yr,color=colors[n],lw=1.6,label=r"Figure, $\eta=0.013$")
        if status == "reference only":
            print(f"{n}: figure q=32 {interp_line(xr,yr,32):.6f}; reference only")
            continue
        y = np.array([float(r["aegis_compensated"]) for r in group])
        xg, yg = curve(q,y)
        rr = max(abs(y/ref-1))
        slope = (np.log(interp_line(xg,yg,16) / interp_line(xg,yg,4)) / np.log(4) - 5/3)
        r_slope = (np.log(interp_line(xr,yr,16) / interp_line(xr,yr,4)) / np.log(4) - 5/3)
        print(f"{source:<37} {status:13} max fractional difference={rr:.4g} "
              f"q=32 Aegis/figure={interp_line(xg,yg,32):.6f}/{interp_line(xr,yr,32):.6f}; "
              f"M(k) slope q=4..16 Aegis/figure={slope:.4f}/{r_slope:.4f}")
        if status == "matched":
            ax.scatter(q,y,s=20,facecolors="white",edgecolors=colors[n],linewidths=1.05,zorder=5)
            if n == 128:
                small.scatter(q,y,s=20,facecolors="white",edgecolors=colors[n],zorder=5,
                              label=r"Aegis, $\eta=0.013$")
        else:
            small.plot(xg,yg,color="#9a456b",ls="--",lw=1.5)
            small.scatter(q,y,s=21,marker="s",facecolors="white",edgecolors="#9a456b",zorder=5,
                          label=rf"Aegis, $\eta={float(group[0]['aegis_eta']):g}$")
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([],[],color="#2a3239",marker="o",markerfacecolor="white",ls="none"))
    labels.append("Aegis shells")
    ax.legend(handles, labels, loc="lower left", ncols=2,fontsize=8,frameon=False,
              columnspacing=0.7,handlelength=2)
    small.legend(loc="lower left",frameon=False,fontsize=8.5)
    fig.savefig(output, metadata={"Title":"Last-snapshot magnetic energy spectra in Aegis"})
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--cache",type=Path,help="Original five highres_magnetic_spectra .h5 caches")
    source.add_argument("--table",type=Path,help="Replot the exported per-shell CSV")
    p.add_argument("--measurements",type=Path,help="Directory written by measure.jl")
    p.add_argument("--out",type=Path,required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.table:
        with args.table.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
    else:
        if not args.measurements:
            p.error("--cache requires --measurements")
        rows = assemble(args.cache, args.measurements)
        with (args.out/"magnetic_shells.csv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]),lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    render(rows, args.out/"magnetic_spectra_overlay.pdf")


if __name__ == "__main__":
    main()
