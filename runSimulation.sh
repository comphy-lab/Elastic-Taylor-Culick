#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/params.sh
source "$SCRIPT_DIR/scripts/params.sh"

usage() {
  cat <<'EOF'
Usage: bash runSimulation.sh [OPTIONS]

Compile and run one Taylor--Culick case.

Options:
  --case FILE       Case source (default: simulationCases/TaylorCulick.c)
  --input FILE      Parameter file (default: default.params)
  --compile-only    Compile the case but do not execute it
  --dry-run         Print the case contract without changing files
  -h, --help        Show this help
EOF
}

resolve_path() {
  local candidate="$1"
  if [[ -f "$candidate" ]]; then
    printf '%s\n' "$candidate"
  elif [[ -f "$SCRIPT_DIR/$candidate" ]]; then
    printf '%s\n' "$SCRIPT_DIR/$candidate"
  else
    return 1
  fi
}

CASE_FILE="simulationCases/TaylorCulick.c"
PARAM_FILE="default.params"
COMPILE_ONLY=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case)
      [[ $# -ge 2 ]] || { echo "ERROR: --case needs a file." >&2; exit 1; }
      CASE_FILE="$2"
      shift 2
      ;;
    --input)
      [[ $# -ge 2 ]] || { echo "ERROR: --input needs a file." >&2; exit 1; }
      PARAM_FILE="$2"
      shift 2
      ;;
    --compile-only) COMPILE_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

CASE_SOURCE="$(resolve_path "$CASE_FILE")" || {
  echo "ERROR: case source not found: $CASE_FILE" >&2
  exit 1
}
PARAM_SOURCE="$(resolve_path "$PARAM_FILE")" || {
  echo "ERROR: parameter file not found: $PARAM_FILE" >&2
  exit 1
}

CASE_NO="$(get_param_value CaseNo "$PARAM_SOURCE")"
if [[ ! "$CASE_NO" =~ ^[0-9]+$ || "$CASE_NO" -lt 1000 ]]; then
  echo "ERROR: CaseNo must be an integer >= 1000: $CASE_NO" >&2
  exit 1
fi

CASE_DIR="$SCRIPT_DIR/simulationCases/c$CASE_NO"
CASE_NAME="$(basename "$CASE_SOURCE" .c)"
LOCAL_SOURCE="$CASE_DIR/$CASE_NAME.c"
EXECUTABLE="$CASE_DIR/$CASE_NAME"
LOCAL_PARAMS="$CASE_DIR/case.params"

echo "CaseNo:       $CASE_NO"
echo "Source:       $CASE_SOURCE"
echo "Parameters:   $PARAM_SOURCE"
echo "Output:       $CASE_DIR"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry run:       mkdir/copy/compile/run skipped."
  exit 0
fi

QCC_BIN="$(command -v qcc || true)"
if [[ -z "$QCC_BIN" ]]; then
  echo "ERROR: qcc not found in PATH." >&2
  exit 1
fi

mkdir -p "$CASE_DIR"
cp "$CASE_SOURCE" "$LOCAL_SOURCE"
cp "$PARAM_SOURCE" "$LOCAL_PARAMS"

echo "Compiling with: $QCC_BIN -O2 -Wall -disable-dimensions"
(
  cd "$CASE_DIR"
  "$QCC_BIN" -O2 -Wall -disable-dimensions \
    -I"$SCRIPT_DIR/src-local" "$CASE_NAME.c" -o "$CASE_NAME" -lm
)

if [[ "$COMPILE_ONLY" -eq 1 ]]; then
  echo "Compilation completed; execution skipped."
  exit 0
fi

echo "Running: $EXECUTABLE case.params"
(
  cd "$CASE_DIR"
  exec "./$CASE_NAME" case.params
)
