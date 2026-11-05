#!/bin/bash
set -eu
SCRIPT_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
nohup bash "$SCRIPT_HOME/run_autonomous.sh" "$@" >> "$PWD/wannier_autonomous.log" 2>&1 < /dev/null &
echo "Wannier runner started (PID $!). Log: $PWD/wannier_autonomous.log"
