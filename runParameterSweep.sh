#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
# shellcheck source=scripts/params.sh
source "$SCRIPT_DIR/scripts/params.sh"

usage() {
  cat <<'EOF'
Usage: bash runParameterSweep.sh [PARAM_FILE] [--dry-run]

Generate the Cartesian product of SWEEP_* entries and run each case through
runSimulation.sh. CASE_START and CASE_END must match the combination count.

Options:
  --dry-run         Print generated cases without compiling or running them
  -h, --help        Show this help
EOF
}

PARAM_FILE="sweep.params"
PARAM_SET=0
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    --)
      shift
      break
      ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      [[ "$PARAM_SET" -eq 0 ]] || {
        echo "ERROR: only one parameter file is accepted." >&2
        exit 1
      }
      PARAM_FILE="$1"
      PARAM_SET=1
      shift
      ;;
  esac
done
[[ "$#" -eq 0 ]] || {
  echo "ERROR: unexpected arguments: $*" >&2
  exit 1
}

if [[ ! -f "$PARAM_FILE" && -f "$SCRIPT_DIR/$PARAM_FILE" ]]; then
  PARAM_FILE="$SCRIPT_DIR/$PARAM_FILE"
fi
[[ -f "$PARAM_FILE" ]] || {
  echo "ERROR: file not found: $PARAM_FILE" >&2
  exit 1
}

CASE_START="$(get_param_value CASE_START "$PARAM_FILE")"
CASE_END="$(get_param_value CASE_END "$PARAM_FILE")"
if [[ ! "$CASE_START" =~ ^[0-9]+$ ||
      ! "$CASE_END" =~ ^[0-9]+$ ||
      "$CASE_START" -lt 1000 ||
      "$CASE_END" -lt "$CASE_START" ]]; then
  echo "ERROR: CASE_START/CASE_END must be integers with 1000 <= start <= end." >&2
  exit 1
fi

declare -a SWEEP_KEYS=()
declare -a SWEEP_SPECS=()
while IFS= read -r key; do
  [[ "$key" == SWEEP_* ]] || continue
  spec="$(get_param_value "$key" "$PARAM_FILE")"
  [[ -n "$spec" ]] || {
    echo "ERROR: empty sweep entry: $key" >&2
    exit 1
  }
  SWEEP_KEYS+=("${key#SWEEP_}")
  SWEEP_SPECS+=("$spec")
done < <(get_param_keys "$PARAM_FILE")

[[ "${#SWEEP_KEYS[@]}" -gt 0 ]] || {
  echo "ERROR: no SWEEP_* entries found in $PARAM_FILE" >&2
  exit 1
}

COMBINATION_COUNT=1
for spec in "${SWEEP_SPECS[@]}"; do
  IFS=',' read -r -a values <<< "$spec"
  [[ "${#values[@]}" -gt 0 ]] || {
    echo "ERROR: empty sweep values." >&2
    exit 1
  }
  COMBINATION_COUNT=$((COMBINATION_COUNT * ${#values[@]}))
done

EXPECTED_COUNT=$((CASE_END - CASE_START + 1))
if [[ "$COMBINATION_COUNT" -ne "$EXPECTED_COUNT" ]]; then
  echo "ERROR: CASE_START..CASE_END describes $EXPECTED_COUNT cases,"
  echo "       but SWEEP_* entries generate $COMBINATION_COUNT combinations." >&2
  exit 1
fi

echo "Sweep file:        $PARAM_FILE"
echo "Combination count: $COMBINATION_COUNT"
echo "Case range:        $CASE_START..$CASE_END"

TEMP_ROOT="$(mktemp -d "/tmp/elastic-taylor-culick-sweep.XXXXXX")"
trap 'rm -rf "$TEMP_ROOT"' EXIT

generate_cases() {
  local index="$1"
  local base_file="$2"
  local temporary_params
  local key
  local spec
  local value
  local case_no
  local -a values

  if [[ "$index" -eq "${#SWEEP_KEYS[@]}" ]]; then
    case_no="$NEXT_CASE"
    NEXT_CASE=$((NEXT_CASE + 1))
    temporary_params="$TEMP_ROOT/case-$case_no.params"
    cp "$base_file" "$temporary_params"
    set_param_in_file CaseNo "$case_no" "$temporary_params"

    echo "Case $case_no:"
    for ((j = 0; j < ${#SWEEP_KEYS[@]}; j++)); do
      key="${SWEEP_KEYS[$j]}"
      echo "  $key=$(get_param_value "$key" "$temporary_params")"
    done

    if [[ "$DRY_RUN" -eq 0 ]]; then
      bash "$SCRIPT_DIR/runSimulation.sh" --input "$temporary_params"
    fi
    return
  fi

  key="${SWEEP_KEYS[$index]}"
  spec="${SWEEP_SPECS[$index]}"
  IFS=',' read -r -a values <<< "$spec"
  for value in "${values[@]}"; do
    value="$(trim_value "$value")"
    temporary_params="$TEMP_ROOT/partial-$index-$NEXT_CASE.params"
    cp "$base_file" "$temporary_params"
    set_param_in_file "$key" "$value" "$temporary_params"
    generate_cases $((index + 1)) "$temporary_params"
  done
}

BASE_PARAMS="$TEMP_ROOT/base.params"
copy_without_sweep_params "$PARAM_FILE" "$BASE_PARAMS"
NEXT_CASE="$CASE_START"
generate_cases 0 "$BASE_PARAMS"
