#!/bin/sh
set -eu

# The only shell entry point; all operations live in the Rust tool.
HORB_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$HORB_ROOT"
exec cargo run --quiet --locked --manifest-path "$HORB_ROOT/Cargo.toml" -p horb-tool -- "$@"
