#!/bin/bash
# Compatibility alias; use ../run_autonomous.sh for new runs.
set -eu
SCRIPT_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$SCRIPT_HOME/run_autonomous.sh" "$@"
