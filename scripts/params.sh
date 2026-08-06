#!/usr/bin/env bash

# Shared shell helpers for CoMPhy key=value parameter files.

get_param_value() {
  local key="$1"
  local file="$2"

  awk -F '=' -v wanted_key="$key" '
    /^[[:space:]]*#/ { next }
    {
      current_key = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", current_key)
      if (current_key != wanted_key)
        next
      value = substr($0, index($0, "=") + 1)
      sub(/[[:space:]]*#.*/, "", value)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      print value
      exit
    }
  ' "$file"
}

get_param_keys() {
  local file="$1"

  awk -F '=' '
    /^[[:space:]]*#/ { next }
    {
      key = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", key)
      if (key != "")
        print key
    }
  ' "$file"
}

set_param_in_file() {
  local key="$1"
  local value="$2"
  local file="$3"
  local temporary_file

  temporary_file="$(mktemp "${TMPDIR:-/tmp}/elastic-taylor-culick-param.XXXXXX")"
  awk -v wanted_key="$key" -v new_value="$value" '
    BEGIN { updated = 0 }
    /^[[:space:]]*#/ { print; next }
    {
      current_key = $1
      sub(/[[:space:]]*=.*/, "", current_key)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", current_key)
      if (current_key == wanted_key) {
        print wanted_key " = " new_value
        updated = 1
      }
      else
        print
    }
    END {
      if (!updated)
        print wanted_key " = " new_value
    }
  ' "$file" > "$temporary_file"
  mv "$temporary_file" "$file"
}

copy_without_sweep_params() {
  local source_file="$1"
  local destination_file="$2"

  awk -F '=' '
    /^[[:space:]]*#/ { print; next }
    {
      key = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", key)
      if (key == "CASE_START" || key == "CASE_END" ||
          key ~ /^SWEEP_/) next
      print
    }
  ' "$source_file" > "$destination_file"
}

trim_value() {
  printf '%s\n' "$1" | awk '{ gsub(/^[[:space:]]+|[[:space:]]+$/, ""); print }'
}
