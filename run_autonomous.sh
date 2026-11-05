#!/bin/bash
set -eu
SCRIPT_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${QE_ENV_FILE:-}" ]; then source "$QE_ENV_FILE"; fi
exec python3 "$SCRIPT_HOME/workflow/wannier_workflow.py" "$@"
