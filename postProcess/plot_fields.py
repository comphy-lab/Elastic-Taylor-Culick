#!/usr/bin/env python3
"""Render one PNG per snapshot from run_postprocess.sh's output and
assemble them into an mp4 with ffmpeg.

Two panels stacked per frame, both mirrored about the midplane y=0 (the
same field on both halves -- the simulated half-domain is symmetric, so
mirroring shows the physical full sheet cross-section):

  wide       velocity magnitude (`Blues`, white at |u|=0 -- the
             `DropsAtLubis`/`Basilisk-101` convention) over the fixed
             `[0, wide-fraction*Ldomain]` window `run_postprocess.sh` wrote.
  comoving   viscous dissipation rate per unit volume (`hot_r`, `LogNorm`
             on the raw value -- `get_fields` writes raw PHI for exactly
             this reason) over the window centred on the tip.

Styling (colormaps, LogNorm dissipation, magenta interface, horizontal
top colorbar) matches
`Taylor-Culick-FEM/postProcess/animate_planar.py`'s rendered style.

Reads ``<case-dir>/postprocess/{facets,fields,wide_fields}/snapshot-*.dat``
(written by ``get_facets``/``get_fields``) and ``case.params`` (for ``h0``).
Requires matplotlib and ffmpeg on PATH -- the ``elastic-tc-postprocess``
conda env has both.

Real LaTeX (``text.usetex``) is the default now that Stokes has a working
system `texlive-latex-extra` install (conda-forge's own `texlive-core` is
missing `latex.ltx` entirely -- do not reintroduce it into this conda env,
which would put a broken `latex` back ahead of `/usr/bin/latex` on PATH).
Static text (axis/colorbar labels) is cached across frames by matplotlib's
TexManager; only the per-frame time label recompiles, so this stays fast
across a full run. Pass --mathtext to force the no-subprocess renderer
instead (needed if this is ever parallelised with multiprocessing, per the
documented LaTeX-daemon deadlock).

Usage
-----
    python3 postProcess/plot_fields.py --case-dir <dir> [--viscoelastic]
                                        [--fps 30] [--mathtext]
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "serif"


def configure_typography(mathtext: bool) -> None:
    if mathtext:
        matplotlib.rcParams["mathtext.fontset"] = "cm"
    else:
        matplotlib.rcParams["font.serif"] = ["Computer Modern Roman"]
        matplotlib.rcParams["text.usetex"] = True
        matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"


import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import numpy as np  # noqa: E402

V_TC = np.sqrt(2.0)
INTERFACE_COLOR = "#FF00C8"


def read_facets(fname: str) -> list[np.ndarray]:
    """Return a list of (2,2) arrays, one per interface segment.

    ``output_facets()`` writes each segment as two `x y` lines followed by
    a blank line; a naive single ``plt.plot`` across the whole file would
    wrongly draw a connecting line between unrelated segments, so segments
    are kept separate and drawn with a `LineCollection`.
    """
    segments = []
    pts: list[tuple[float, float]] = []
    with open(fname) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                if len(pts) == 2:
                    segments.append(np.array(pts))
                pts = []
                continue
            x, y = map(float, line.split())
            pts.append((x, y))
    if len(pts) == 2:
        segments.append(np.array(pts))
    return segments


def read_fields(fname: str):
    """Return columns as 2D arrays shaped (nx, ny) from the
    ``# nx N ny N`` header `get_fields` writes."""
    with open(fname) as fh:
        header = fh.readline()
    m = re.match(r"#\s*nx\s+(\d+)\s+ny\s+(\d+)", header)
    if not m:
        raise ValueError(f"{fname}: missing '# nx .. ny ..' header")
    nx, ny = int(m.group(1)), int(m.group(2))
    data = np.loadtxt(fname, comments="#")
    ncols = data.shape[1]
    names = ["x", "y", "f", "ux", "uy", "kappa", "phi"]
    if ncols > len(names):
        names += ["A11", "A12", "A22", "T11", "T12", "T22"][: ncols - len(names)]
    cols = {name: data[:, i].reshape(nx, ny) for i, name in enumerate(names)}
    # interpolate() returns Basilisk's `nodata` (~1e30) for a point outside
    # the restored domain -- run_postprocess.sh clamps its windows to avoid
    # this, but mask it here too so a future out-of-domain window degrades
    # to a visibly-blank patch instead of a wrecked colour range (a single
    # 1e30 dominates a 99.5th-percentile or a LogNorm span).
    for name in ("ux", "uy", "kappa", "phi"):
        cols[name] = np.where(np.abs(cols[name]) > 1e29, np.nan, cols[name])
    return cols


def phi_range(phi: np.ndarray) -> tuple[float, float]:
    # Same construction as animate_planar.py: high end from the 99.5th
    # percentile (robust to a few near-singular interface cells), low end
    # four decades below so the colormap still resolves the far-field decay.
    # nanpercentile: out-of-domain cells are masked to NaN in read_fields().
    hi = max(float(np.nanpercentile(phi, 99.5)), 1e-6)
    return (max(hi * 1e-4, 1e-8), hi)


def mirrored_panel(ax, cax, cols, segments, values, cmap, norm_kwargs,
                    cbar_label, time_label=None):
    """Plot ``values`` mirrored about y=0 on both halves of ``ax``, with a
    single horizontal colorbar on ``cax`` above it."""
    mesh = ax.pcolormesh(cols["x"], cols["y"], values, cmap=cmap,
                          shading="gouraud", **norm_kwargs)
    ax.pcolormesh(cols["x"], -cols["y"], values, cmap=cmap,
                   shading="gouraud", **norm_kwargs)

    mirrored = [np.column_stack([s[:, 0], -s[:, 1]]) for s in segments]
    ax.add_collection(LineCollection(segments + mirrored,
                                      colors=INTERFACE_COLOR, linewidths=1.6,
                                      zorder=5, capstyle="round"))
    ax.axhline(0.0, color="0.45", lw=1.0, ls=(0, (6, 4)), zorder=6)

    ax.set_xlim(cols["x"].min(), cols["x"].max())
    ax.set_ylim(-cols["y"].max(), cols["y"].max())
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/h_0$", fontsize=24, labelpad=7)
    ax.set_ylabel(r"$y/h_0$", fontsize=24, labelpad=7)
    ax.tick_params(which="both", direction="out", width=2, labelsize=18,
                    pad=5)
    for sp in ax.spines.values():
        sp.set_linewidth(2)

    if time_label is not None:
        box = dict(boxstyle="round,pad=0.28", fc="k", ec="none", alpha=0.55)
        ax.text(0.985, 0.90, time_label, transform=ax.transAxes, ha="right",
                 fontsize=22, color="w", bbox=box, zorder=8)

    cb = plt.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.set_label(cbar_label, fontsize=22, labelpad=10)
    cb.ax.tick_params(labelsize=17, pad=4)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.xaxis.set_ticks_position("top")
    return mesh


def render_frame(facets_file: str, wide_file: str, comoving_file: str,
                  out_png: str, t: float) -> None:
    segments = read_facets(facets_file)
    wide = read_fields(wide_file)
    com = read_fields(comoving_file)
    umag_wide = np.hypot(wide["ux"], wide["uy"]) / V_TC
    phi_lim = phi_range(com["phi"])

    # ---- geometry: panel heights derived from each row's (mirrored) data
    # aspect at a shared panel width -- a fixed height ratio plus
    # set_aspect("equal") leaves the differently-shaped panel padded with
    # blank canvas rather than filling its box.
    dx_wide = wide["x"].max() - wide["x"].min()
    dy_wide = wide["y"].max() - wide["y"].min()  # half-span; mirrored -> 2x
    dx_com = com["x"].max() - com["x"].min()
    dy_com = com["y"].max() - com["y"].min()

    PW = 9.0                                   # shared panel width, inches
    RH1 = PW * (2 * dy_wide) / dx_wide          # row 1 height (mirrored)
    RH2 = PW * (2 * dy_com) / dx_com            # row 2 height (mirrored)

    ml, mr = 1.30, 0.30
    cb_h = 0.30
    row_gap = 1.55                              # clears row 1's x tick/label
    gap_hi = 0.32                               # colorbar sits just above
    mb, mt = 1.10, 0.30

    FW = ml + PW + mr
    FH = mt + cb_h + gap_hi + RH1 + row_gap + cb_h + gap_hi + RH2 + mb

    fig = plt.figure(figsize=(FW, FH))

    y_row2 = mb / FH
    y_row2_cb = (mb + RH2 + gap_hi) / FH
    y_row1 = (mb + RH2 + gap_hi + cb_h + row_gap) / FH
    y_row1_cb = (mb + RH2 + gap_hi + cb_h + row_gap + RH1 + gap_hi) / FH

    ax_wide = fig.add_axes([ml / FW, y_row1, PW / FW, RH1 / FH])
    cax_wide = fig.add_axes([ml / FW, y_row1_cb, PW / FW, cb_h / FH])
    ax_com = fig.add_axes([ml / FW, y_row2, PW / FW, RH2 / FH])
    cax_com = fig.add_axes([ml / FW, y_row2_cb, PW / FW, cb_h / FH])

    t_label = rf"$t = {t:5.2f}$"

    mirrored_panel(ax_wide, cax_wide, wide, segments, umag_wide, "Blues",
                    dict(vmin=0.0, vmax=1.0),
                    r"$|u|\,/\,V_{\mathrm{TC}}$", time_label=t_label)

    # Mark row 2's comoving window on row 1, so the two panels' relationship
    # stays legible even when the rim is a small feature of the wide view.
    com_xmin, com_xmax = com["x"].min(), com["x"].max()
    com_ymax = com["y"].max()
    ax_wide.add_patch(Rectangle(
        (com_xmin, -com_ymax), com_xmax - com_xmin, 2 * com_ymax,
        facecolor="none", edgecolor="0.4", linewidth=1.5,
        linestyle=(0, (4, 3)), zorder=7))

    mirrored_panel(ax_com, cax_com, com, segments, com["phi"], "hot_r",
                    dict(norm=LogNorm(*phi_lim)),
                    r"$\varepsilon\,/\,(\sigma/h_0)\sqrt{\sigma/\rho h_0}$",
                    time_label=t_label)

    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--viscoelastic", action="store_true",
                     help="unused here beyond validating the fields file "
                          "has the expected extra columns; kept for "
                          "symmetry with run_postprocess.sh's flag")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--mathtext", action="store_true",
                     help="use matplotlib's mathtext instead of real LaTeX "
                          "-- only needed if this is ever run under "
                          "multiprocessing (LaTeX-daemon deadlock)")
    args = ap.parse_args()
    configure_typography(args.mathtext)

    case_dir = args.case_dir
    pp_dir = os.path.join(case_dir, "postprocess")
    facets_dir = os.path.join(pp_dir, "facets")
    fields_dir = os.path.join(pp_dir, "fields")
    wide_dir = os.path.join(pp_dir, "wide_fields")
    frames_dir = os.path.join(pp_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    fields_files = sorted(glob.glob(os.path.join(fields_dir, "snapshot-*.dat")))
    if not fields_files:
        raise SystemExit(f"no snapshot-*.dat under {fields_dir}; "
                          f"run run_postprocess.sh first")

    for fields_file in fields_files:
        tstamp = os.path.basename(fields_file)[len("snapshot-"):-len(".dat")]
        facets_file = os.path.join(facets_dir, f"snapshot-{tstamp}.dat")
        wide_file = os.path.join(wide_dir, f"snapshot-{tstamp}.dat")
        out_png = os.path.join(frames_dir, f"frame-{tstamp}.png")
        render_frame(facets_file, wide_file, fields_file, out_png,
                     float(tstamp))
        print(f"  {out_png}")

    video_path = os.path.join(pp_dir, "taylor_culick.mp4")
    cmd = [
        "ffmpeg", "-y", "-framerate", str(args.fps),
        "-pattern_type", "glob", "-i", os.path.join(frames_dir, "frame-*.png"),
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v", "libx264", "-pix_fmt", "yuv420p", video_path,
    ]
    subprocess.run(cmd, check=True)
    print(f"-> {video_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
