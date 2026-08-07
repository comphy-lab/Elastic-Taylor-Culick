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
  --outdir DIR      Run directory (default: simulationCases/c<CaseNo>).
                    Use this to keep run data out of the checkout.
  --openmp          Compile with -fopenmp (honours OMP_NUM_THREADS)
  --mpi             Compile for MPI (-D_MPI=1 -D_DEFAULT_SOURCE)
  --np N            Launch the MPI build under `mpirun -n N` (needs --mpi)
  --compile-only    Compile the case but do not execute it
  --dry-run         Print the case contract without changing files
  -h, --help        Show this help

Build a parallel binary.  Basilisk arms floating-point exception trapping
in grid/config.h under `(_GNU_SOURCE || __APPLE__) && !_OPENMP`, and qcc
passes -D_GNU_SOURCE itself, so a plain serial build of these cases has
the trap armed and dies with SIGFPE on the first timestep.  An OpenMP
build defines _OPENMP and compiles the trap out; for MPI, qcc strips
-D_GNU_SOURCE from $CC99 (qcc.c:415), which is why --mpi passes
-D_DEFAULT_SOURCE -- it still exposes MAP_ANONYMOUS and MADV_DONTNEED but
leaves enable_fpe() a no-op, matching the OpenMP build.  Never pass
-D_GNU_SOURCE to an MPI build; qcc cannot combine MPI with OpenMP, so
that build has no _OPENMP to disarm the trap.
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
OUT_DIR=""
OPENMP=0
USE_MPI=0
NP=""
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
    --outdir)
      [[ $# -ge 2 ]] || { echo "ERROR: --outdir needs a directory." >&2; exit 1; }
      OUT_DIR="$2"
      shift 2
      ;;
    --openmp) OPENMP=1; shift ;;
    --mpi) USE_MPI=1; shift ;;
    --np)
      [[ $# -ge 2 ]] || { echo "ERROR: --np needs a rank count." >&2; exit 1; }
      NP="$2"
      shift 2
      ;;
    --compile-only) COMPILE_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ "$USE_MPI" -eq 1 && "$OPENMP" -eq 1 ]]; then
  echo "ERROR: --mpi and --openmp are mutually exclusive; qcc reports" >&2
  echo "       'OpenMP cannot be used with MPI (yet)' and drops -fopenmp." >&2
  exit 1
fi
if [[ -n "$NP" ]]; then
  if [[ ! "$NP" =~ ^[0-9]+$ || "$NP" -lt 1 ]]; then
    echo "ERROR: --np must be a positive integer: $NP" >&2
    exit 1
  fi
  if [[ "$USE_MPI" -eq 0 ]]; then
    echo "ERROR: --np requires --mpi." >&2
    exit 1
  fi
fi

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

if [[ -n "$OUT_DIR" ]]; then
  mkdir -p "$OUT_DIR"
  CASE_DIR="$(cd "$OUT_DIR" && pwd)"
else
  CASE_DIR="$SCRIPT_DIR/simulationCases/c$CASE_NO"
fi
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

QCC_FLAGS=(-O2 -Wall -disable-dimensions)
if [[ "$OPENMP" -eq 1 ]]; then
  QCC_FLAGS+=(-fopenmp)
fi
if [[ "$USE_MPI" -eq 1 ]]; then
  # -D_DEFAULT_SOURCE, never -D_GNU_SOURCE: see --help.
  QCC_FLAGS+=(-D_MPI=1 -D_DEFAULT_SOURCE)
  # qcc reads its backend compiler from $CC99 (qcc.c:136); the default one
  # does not know where mpi.h lives, so an unset CC99 fails at grid/config.h.
  if [[ -z "${CC99:-}" ]]; then
    if ! command -v mpicc >/dev/null 2>&1; then
      echo "ERROR: --mpi needs mpicc in PATH, or CC99 set explicitly." >&2
      exit 1
    fi
    export CC99="mpicc -std=c99"
  fi
  echo "MPI backend:  CC99=$CC99"
fi

if [[ "$OPENMP" -eq 0 && "$USE_MPI" -eq 0 && "$COMPILE_ONLY" -eq 0 ]]; then
  echo "WARNING: building serial (neither --openmp nor --mpi)." >&2
  echo "         qcc passes -D_GNU_SOURCE and there is no _OPENMP to disarm" >&2
  echo "         Basilisk's FP-exception trap, so these cases abort with" >&2
  echo "         SIGFPE on the first timestep.  Use --openmp or --mpi." >&2
fi

echo "Compiling with: $QCC_BIN ${QCC_FLAGS[*]}"
(
  cd "$CASE_DIR"
  "$QCC_BIN" "${QCC_FLAGS[@]}" \
    -I"$SCRIPT_DIR/src-local" "$CASE_NAME.c" -o "$CASE_NAME" -lm
)

if [[ "$COMPILE_ONLY" -eq 1 ]]; then
  echo "Compilation completed; execution skipped."
  exit 0
fi

if [[ "$USE_MPI" -eq 1 && -n "$NP" ]]; then
  if ! command -v mpirun >/dev/null 2>&1; then
    echo "ERROR: --np given but mpirun is not in PATH." >&2
    exit 1
  fi
  echo "Running: mpirun -n $NP $EXECUTABLE case.params"
  cd "$CASE_DIR"
  exec mpirun -n "$NP" "./$CASE_NAME" case.params
fi

echo "Running: $EXECUTABLE case.params"
(
  cd "$CASE_DIR"
  exec "./$CASE_NAME" case.params
)
