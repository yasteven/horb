#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-cat_all_scripts.txt}"
DIR="${2:-scripts}"

if [[ ! -d "$DIR" ]]; then
  echo "missing directory: $DIR" >&2
  exit 1
fi

{
  echo "# cat_all_scripts.txt"
  echo "# generated: $(date -Is)"
  echo "# repo: $(pwd)"
  echo "# dir: $DIR"
  echo

  find "$DIR" -maxdepth 1 -type f | sort | while read -r f; do
    echo
    echo "================================================================================"
    echo "FILE: $f"
    echo "SIZE: $(wc -c < "$f") bytes"
    echo "LINES: $(wc -l < "$f")"
    echo "MODE: $(stat -c '%A %a' "$f")"
    echo "================================================================================"
    echo

    case "$f" in
      *.sh|*.py)
        cat "$f"
        ;;
      *)
        echo "[skipped: not .sh or .py]"
        ;;
    esac

    echo
  done
} > "$OUT"

echo "wrote $OUT"
ls -lh "$OUT"
