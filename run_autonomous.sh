#!/bin/bash
set -eu
SCRIPT_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$SCRIPT_HOME/run_wannier_autonomous.sh" "$@"
