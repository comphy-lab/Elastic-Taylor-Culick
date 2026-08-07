#!/usr/bin/env bash
#
# Post-process one run directory's dump snapshots directly: compile the
# Basilisk snapshot readers in this folder, then for every
# intermediate/snapshot-* file extract interface facets, an independent
# tip-position/velocity cross-check, and two field grids: one centred on
# the retracting tip (moves with it every snapshot) and one over a fixed
# [0, wide-fraction*Ldomain] span of the domain.
#
# The comoving window is always an exact square centred on the tip, with
# half-width scaled off the rim -- see the R_est derivation below. A fixed
# 1:1 data aspect every snapshot is load-bearing, not cosmetic: plot_fields.py
# derives the comoving panel's PIXEL SIZE from the data aspect ratio, so a
# window whose shape changes frame to frame produces a differently-sized PNG
# every frame -- ffmpeg then either fails or produces a visibly jumpy video.
# Squaring it makes the panel size constant by construction.
#
# Usage:
#   run_postprocess.sh --case-dir DIR [--viscoelastic]
#                       [--window-coef N] [--window-offset N] [--ny N]
#                       [--wide-fraction F] [--wide-y N] [--wide-ny N]
#                       [--video] [--fps N]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../scripts/params.sh
source "$REPO_ROOT/scripts/params.sh"

usage() {
  cat <<'EOF'
Usage: run_postprocess.sh --case-dir DIR [OPTIONS]

Options:
  --case-dir DIR      Run directory containing case.params and
                       intermediate/snapshot-* (required)
  --viscoelastic       Also extract A11/A12/A22/T11/T12/T22 (the snapshot
                       must come from TaylorCulickPlanar.c, not the
                       Newtonian-only case)
  --window-coef N      Comoving window half-width W = coef*R_est + offset,
                       in units of h0, applied symmetrically as
                       [x_tip - W, x_tip + W] x [-W, W] (default coef 1.8)
  --window-offset N    (default offset 2.0) -- R_est =
                       sqrt(h0*(x_tip - xtip0)/pi) is a mass-conservation
                       estimate of the rim's radius: the liquid swept by
                       the retracting edge (h0*(x_tip - xtip0) per unit
                       span, upper half only) is assumed to collect into a
                       half-disk of radius R_est. Checked against a real
                       run's measured rim apex height at t=10/20/30/40:
                       tracks R_est to within 11% (3% by t=40). Coef/offset
                       are tuned so the blob fills a substantial, roughly
                       constant fraction of the square at every time (a
                       wider window showed more of the neck-to-film
                       reflattening but made the blob look lost in mostly-
                       white canvas, worsened by hot_r's low end being
                       near-white with no axis border to give scale).
                       R_est needs no smoothing (it is a monotonic function
                       of x_tip, already smooth) and degrades gracefully at
                       t=0 (R_est=0, so the offset alone sets the window),
                       so no separate floor logic is needed. Not clamped to
                       the domain: an unclamped x_tip - W*h0 can go
                       negative early in a run, but interpolate() returns
                       Basilisk's `nodata` (~1e30) for points outside the
                       domain rather than erroring, and plot_fields.py
                       masks that to NaN -- clamping would break the exact
                       symmetry the fixed panel size depends on.
  --ny N               Grid rows in the comoving field window (default 200
                       -- scaled up from the original 80 to hold ~0.1*h0
                       resolution at the larger comoving windows)
  --wide-fraction F    Wide-view window as a fraction of Ldomain, fixed at
                       [0, F*Ldomain] for every snapshot (default 0.5 -- the
                       point of row 1 is the whole domain's air flow over
                       the run, not a close-up of the rim; row 2's comoving
                       window is marked on row 1 with a grey box so the two
                       panels' relationship stays legible even when the rim
                       is a small feature of row 1)
  --wide-y N           Wide-view window height above the midplane, in
                       units of h0 (default 4 -- the rim apex reaches
                       ~3.9*h0 by the end of a real run, so 3 clipped it;
                       row 1 only needs to show part of the blob, not
                       resolve it, so no need to match row 2's scale)
  --wide-ny N          Grid rows in the wide-view window (default 200)
  --video              Render PNG frames with plot_fields.py and assemble
                       an mp4 with ffmpeg (needs the elastic-tc-postprocess
                       conda env: matplotlib + ffmpeg)
  --fps N              Video frame rate (default 30)
  -h, --help           Show this help
EOF
}

CASE_DIR=""
VISCOELASTIC=0
WINDOW_COEF=1.8
WINDOW_OFFSET=2.0
NY=200
WIDE_FRACTION=0.5
WIDE_Y=4
WIDE_NY=200
VIDEO=0
FPS=30

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-dir) CASE_DIR="$2"; shift 2 ;;
    --viscoelastic) VISCOELASTIC=1; shift ;;
    --window-coef) WINDOW_COEF="$2"; shift 2 ;;
    --window-offset) WINDOW_OFFSET="$2"; shift 2 ;;
    --ny) NY="$2"; shift 2 ;;
    --wide-fraction) WIDE_FRACTION="$2"; shift 2 ;;
    --wide-y) WIDE_Y="$2"; shift 2 ;;
    --wide-ny) WIDE_NY="$2"; shift 2 ;;
    --video) VIDEO=1; shift ;;
    --fps) FPS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$CASE_DIR" ]] || { echo "ERROR: --case-dir is required." >&2; exit 1; }
[[ -d "$CASE_DIR" ]] || { echo "ERROR: no such directory: $CASE_DIR" >&2; exit 1; }
[[ -f "$CASE_DIR/case.params" ]] || {
  echo "ERROR: $CASE_DIR/case.params not found." >&2; exit 1; }
SNAP_DIR="$CASE_DIR/intermediate"
[[ -d "$SNAP_DIR" ]] || {
  echo "ERROR: $SNAP_DIR not found (no snapshots to process)." >&2; exit 1; }

QCC="$(which qcc || true)"
[[ -n "$QCC" ]] || { echo "ERROR: qcc not found on PATH." >&2; exit 1; }

H0="$(get_param_value h0 "$CASE_DIR/case.params")"
H0="${H0:-1.0}"
MU1="$(get_param_value mu1 "$CASE_DIR/case.params")"
MU2="$(get_param_value mu2 "$CASE_DIR/case.params")"
[[ -n "$MU1" && -n "$MU2" ]] || {
  echo "ERROR: mu1/mu2 not found in $CASE_DIR/case.params" >&2; exit 1; }
LDOMAIN="$(get_param_value Ldomain "$CASE_DIR/case.params")"
[[ -n "$LDOMAIN" ]] || {
  echo "ERROR: Ldomain not found in $CASE_DIR/case.params" >&2; exit 1; }
XTIP0="$(get_param_value xtip0 "$CASE_DIR/case.params")"
[[ -n "$XTIP0" ]] || {
  echo "ERROR: xtip0 not found in $CASE_DIR/case.params" >&2; exit 1; }
WIDE_XMAX=$(awk -v l="$LDOMAIN" -v f="$WIDE_FRACTION" 'BEGIN{print f*l}')

BUILD_DIR="$SCRIPT_DIR/.build"
mkdir -p "$BUILD_DIR"

echo "Compiling postprocessors with $QCC ..."
# qcc writes a `<source>-cpp.c` intermediate next to the source path as
# given; an absolute source path breaks that when cwd != the source
# directory, so cd here and pass bare filenames.
(
  cd "$SCRIPT_DIR"
  "$QCC" -disable-dimensions -O2 -o "$BUILD_DIR/get_facets" \
    get_facets.c -lm
  "$QCC" -disable-dimensions -O2 -o "$BUILD_DIR/get_tip" \
    get_tip.c -lm
  if [[ "$VISCOELASTIC" -eq 1 ]]; then
    "$QCC" -disable-dimensions -O2 -DVISCOELASTIC=1 \
      -o "$BUILD_DIR/get_fields" get_fields.c -lm
  else
    "$QCC" -disable-dimensions -O2 \
      -o "$BUILD_DIR/get_fields" get_fields.c -lm
  fi
)

OUT_DIR="$CASE_DIR/postprocess"
FACETS_DIR="$OUT_DIR/facets"
FIELDS_DIR="$OUT_DIR/fields"
WIDE_DIR="$OUT_DIR/wide_fields"
mkdir -p "$FACETS_DIR" "$FIELDS_DIR" "$WIDE_DIR"
TIP_CSV="$OUT_DIR/tip_snapshots.csv"
echo "t,x_tip,x_tip_vof,x_tip_global,u_tip_x" > "$TIP_CSV"

# Snapshot filenames are `snapshot-<t>`; sort numerically on the suffix, not
# lexicographically (lexicographic sort would put snapshot-10.0 before
# snapshot-2.0).
mapfile -t SNAPSHOTS < <(
  find "$SNAP_DIR" -maxdepth 1 -name 'snapshot-*' -printf '%f\n' |
    sort -t- -k2 -g
)
[[ ${#SNAPSHOTS[@]} -gt 0 ]] || {
  echo "ERROR: no snapshot-* files under $SNAP_DIR" >&2; exit 1; }

echo "Processing ${#SNAPSHOTS[@]} snapshots from $CASE_DIR ..."
for snap in "${SNAPSHOTS[@]}"; do
  tstamp="${snap#snapshot-}"
  snap_path="$SNAP_DIR/$snap"

  tip_line="$("$BUILD_DIR/get_tip" "$snap_path" "$H0")"
  echo "$tip_line" | tr ' ' ',' >> "$TIP_CSV"
  xtip="$(awk '{print $2}' <<< "$tip_line")"

  "$BUILD_DIR/get_facets" "$snap_path" > "$FACETS_DIR/snapshot-$tstamp.dat"

  # Comoving window: exact square, tip exactly centred -- see the
  # usage-header derivation. R_est is a mass-conservation estimate: liquid
  # swept by the retracting edge collects into a half-disk of this radius.
  # Deliberately NOT clamped to the domain: clamping only one side would
  # break the exact tip-centred symmetry the fixed panel size depends on
  # (see usage header for why interpolate() going out of domain is safe).
  r_est=$(awk -v x="$xtip" -v x0="$XTIP0" -v h="$H0" \
    'BEGIN{d=x-x0; if(d<0)d=0; print sqrt(h*d/3.14159265358979)}')
  w=$(awk -v r="$r_est" -v c="$WINDOW_COEF" -v o="$WINDOW_OFFSET" \
    'BEGIN{print c*r+o}')
  xmin=$(awk -v x="$xtip" -v w="$w" 'BEGIN{print x-w}')
  xmax=$(awk -v x="$xtip" -v w="$w" 'BEGIN{print x+w}')
  "$BUILD_DIR/get_fields" "$snap_path" "$xmin" "$xmax" 0 "$w" "$NY" \
    "$MU1" "$MU2" > "$FIELDS_DIR/snapshot-$tstamp.dat"

  # Wide view is a fixed [0, wide-fraction*Ldomain] window, not tip-tracked.
  wide_ymax=$(awk -v h="$H0" -v w="$WIDE_Y" 'BEGIN{print w*h}')
  "$BUILD_DIR/get_fields" "$snap_path" 0 "$WIDE_XMAX" 0 "$wide_ymax" \
    "$WIDE_NY" "$MU1" "$MU2" > "$WIDE_DIR/snapshot-$tstamp.dat"
done

echo "Wrote $TIP_CSV, $FACETS_DIR/, $FIELDS_DIR/, $WIDE_DIR/"

if [[ "$VIDEO" -eq 1 ]]; then
  echo "Rendering frames and assembling video ..."
  python3 "$SCRIPT_DIR/plot_fields.py" --case-dir "$CASE_DIR" \
    $([[ "$VISCOELASTIC" -eq 1 ]] && echo --viscoelastic) \
    --fps "$FPS"
fi
