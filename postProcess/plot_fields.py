#!/usr/bin/env python3
"""Render one PNG per snapshot from run_postprocess.sh's output and
assemble them into an mp4 with ffmpeg.

Reads ``<case-dir>/postprocess/{facets,fields}/snapshot-*.dat`` (written by
``get_facets``/``get_fields``) and ``case.params`` (for ``h0``, used only in
the title). Requires matplotlib and ffmpeg on PATH -- the
``elastic-tc-postprocess`` conda env has both.

Usage
-----
    python3 postProcess/plot_fields.py --case-dir <dir> [--viscoelastic]
                                        [--fps 30] [--field velocity]
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
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


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
    """Return (x, y, columns) as 2D arrays shaped (nx, ny) from the
    ``# nx N ny N`` header `get_fields` writes, plus the column names."""
    with open(fname) as fh:
        header = fh.readline()
    m = re.match(r"#\s*nx\s+(\d+)\s+ny\s+(\d+)", header)
    if not m:
        raise ValueError(f"{fname}: missing '# nx .. ny ..' header")
    nx, ny = int(m.group(1)), int(m.group(2))
    data = np.loadtxt(fname, comments="#")
    ncols = data.shape[1]
    names = ["x", "y", "f", "ux", "uy", "kappa"]
    if ncols > len(names):
        names += ["A11", "A12", "A22", "T11", "T12", "T22"][: ncols - len(names)]
    cols = {name: data[:, i].reshape(nx, ny) for i, name in enumerate(names)}
    return cols


def render_frame(facets_file: str, fields_file: str, out_png: str,
                  field: str, h0: float) -> None:
    cols = read_fields(fields_file)
    x, y = cols["x"], cols["y"]

    if field == "velocity":
        c = np.hypot(cols["ux"], cols["uy"])
        label = r"$|u|$"
    elif field == "kappa":
        c = np.where(np.abs(cols["kappa"]) < 1e29, cols["kappa"], np.nan)
        label = r"$\kappa$"
    elif field in cols:
        c = cols[field]
        label = field
    else:
        raise ValueError(f"unknown --field {field!r}; have {list(cols)}")

    fig, ax = plt.subplots(figsize=(8, 3), dpi=140)
    mesh = ax.pcolormesh(x, y, np.ma.masked_invalid(c), shading="nearest",
                          cmap="viridis")
    fig.colorbar(mesh, ax=ax, label=label, shrink=0.85)

    segments = read_facets(facets_file)
    if segments:
        ax.add_collection(LineCollection(segments, colors="white",
                                          linewidths=1.2))

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_xlabel("$x / h_0$")
    ax.set_ylabel("$y / h_0$")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def get_param_value(key: str, params_file: str) -> str | None:
    with open(params_file) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            if k.strip() == key:
                return v.strip()
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--viscoelastic", action="store_true",
                     help="unused here beyond validating the fields file "
                          "has the expected extra columns; kept for "
                          "symmetry with run_postprocess.sh's flag")
    ap.add_argument("--field", default="velocity",
                     help="velocity | kappa | f | ux | uy | A11 | ... "
                          "(default: velocity)")
    ap.add_argument("--fps", type=int, default=30)
    args = ap.parse_args()

    case_dir = args.case_dir
    pp_dir = os.path.join(case_dir, "postprocess")
    facets_dir = os.path.join(pp_dir, "facets")
    fields_dir = os.path.join(pp_dir, "fields")
    frames_dir = os.path.join(pp_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    params_file = os.path.join(case_dir, "case.params")
    h0 = float(get_param_value("h0", params_file) or 1.0)

    fields_files = sorted(glob.glob(os.path.join(fields_dir, "snapshot-*.dat")))
    if not fields_files:
        raise SystemExit(f"no snapshot-*.dat under {fields_dir}; "
                          f"run run_postprocess.sh first")

    for fields_file in fields_files:
        tstamp = os.path.basename(fields_file)[len("snapshot-"):-len(".dat")]
        facets_file = os.path.join(facets_dir, f"snapshot-{tstamp}.dat")
        out_png = os.path.join(frames_dir, f"frame-{tstamp}.png")
        render_frame(facets_file, fields_file, out_png, args.field, h0)
        print(f"  {out_png}")

    video_path = os.path.join(pp_dir, f"{args.field}.mp4")
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
