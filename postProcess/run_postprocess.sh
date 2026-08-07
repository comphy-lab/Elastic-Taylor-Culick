#!/usr/bin/env bash
#
# Post-process one run directory's dump snapshots directly: compile the
# Basilisk snapshot readers in this folder, then for every
# intermediate/snapshot-* file extract interface facets, an independent
# tip-position/velocity cross-check, and a field grid windowed around the
# retracting tip (the tip moves O(100) sheet thicknesses over a run, so a
# fixed window would either miss the early interface or crop the late one).
#
# Usage:
#   run_postprocess.sh --case-dir DIR [--viscoelastic]
#                       [--window-behind N] [--window-ahead N] [--window-y N]
#                       [--ny N] [--video] [--fps N]
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
  --window-behind N    Field window extent behind the tip, in units of h0
                       (default 3)
  --window-ahead N     Field window extent ahead of the tip, in units of h0
                       (default 1)
  --window-y N         Field window height above the midplane, in units of
                       h0 (default 2)
  --ny N               Grid rows in the field window (default 80)
  --video              Render PNG frames with plot_fields.py and assemble
                       an mp4 with ffmpeg (needs the elastic-tc-postprocess
                       conda env: matplotlib + ffmpeg)
  --fps N              Video frame rate (default 30)
  -h, --help           Show this help
EOF
}

CASE_DIR=""
VISCOELASTIC=0
WINDOW_BEHIND=3
WINDOW_AHEAD=1
WINDOW_Y=2
NY=80
VIDEO=0
FPS=30

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-dir) CASE_DIR="$2"; shift 2 ;;
    --viscoelastic) VISCOELASTIC=1; shift ;;
    --window-behind) WINDOW_BEHIND="$2"; shift 2 ;;
    --window-ahead) WINDOW_AHEAD="$2"; shift 2 ;;
    --window-y) WINDOW_Y="$2"; shift 2 ;;
    --ny) NY="$2"; shift 2 ;;
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
mkdir -p "$FACETS_DIR" "$FIELDS_DIR"
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

  xmin=$(awk -v x="$xtip" -v h="$H0" -v w="$WINDOW_BEHIND" 'BEGIN{print x-w*h}')
  xmax=$(awk -v x="$xtip" -v h="$H0" -v w="$WINDOW_AHEAD" 'BEGIN{print x+w*h}')
  ymax=$(awk -v h="$H0" -v w="$WINDOW_Y" 'BEGIN{print w*h}')
  "$BUILD_DIR/get_fields" "$snap_path" "$xmin" "$xmax" 0 "$ymax" "$NY" \
    > "$FIELDS_DIR/snapshot-$tstamp.dat"
done

echo "Wrote $TIP_CSV, $FACETS_DIR/, $FIELDS_DIR/"

if [[ "$VIDEO" -eq 1 ]]; then
  echo "Rendering frames and assembling video ..."
  python3 "$SCRIPT_DIR/plot_fields.py" --case-dir "$CASE_DIR" \
    $([[ "$VISCOELASTIC" -eq 1 ]] && echo --viscoelastic) \
    --fps "$FPS"
fi
