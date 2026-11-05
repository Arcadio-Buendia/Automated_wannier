# Running and monitoring

[Documentation index](README.md) · [Project README](../README.md)

## Repository layout and entry point
Use **`./run_autonomous.sh`** for foreground execution or with your own `nohup`/log redirection. It launches the Python workflow directly. The optional `scripts/launch_autonomous.sh` helper starts it in the background; `scripts/run_wannier_autonomous.sh` is an alias for older integrations. The former top-level helper paths have moved into `scripts/`.

```text
run_autonomous.sh          main entry point
workflow/                  Python workflow, analysis, and extension modules
knowledge/                 Jev decision rules and source links
examples/                  fresh, reuse, and parallel example configurations
scripts/                   optional background helper and compatibility alias
tests/                     workflow contract and logging checks
docs/                      detailed guides and reference
README.md
requirements.txt
```

Keep your own `wannier_config.json` and `.env` in the calculation folder. Paths in that config remain relative to its folder; launching from another directory works with `--config /path/to/wannier_config.json`. Repository updates do not move calculation results or scratch data.

## Launch in the foreground or detach
From the calculation folder, with the Python environment active:

```bash
bash ~/tools/Automated_wannier/run_autonomous.sh --config wannier_config.json
```

If the scripts and configuration are in the current folder, `./run_autonomous.sh` is equivalent. Without `--config`, the launcher reads `wannier_config.json` from the current working directory. Use this as the standard entry point. The older name is retained only as `scripts/run_wannier_autonomous.sh`, a compatibility alias.

To keep running after closing the console, start a new job with:

```bash
nohup bash ~/tools/Automated_wannier/run_autonomous.sh \
  --config wannier_config.json > wannier.log 2>&1 < /dev/null &
echo $! > wannier.pid
disown
```

Use `disown` from an interactive Bash shell. Do not start another copy of an already running job. `scripts/launch_autonomous.sh --config ...` is an alternative background launcher; it appends console output to `wannier_autonomous.log` in the current directory. On a scheduled cluster, put the launch command in the site's batch job and use its allocation and launcher.

An analysis-only trial uses:

```bash
bash ~/tools/Automated_wannier/run_autonomous.sh \
  --config wannier_config.json --analyze-only
```

This still copies save data, runs NSCF if required, runs `projwfc`, and produces Python analysis. It makes no Jev requests and performs no localization. Use a **fresh `work_dir` for a later full run**; the program does not resume an analysis-only directory.

## Follow logs and inspect results
For a detached run:

```bash
tail -F wannier.log
```

The first large save copy can take time before stage output appears. The printed scratch path is the copy destination. Discover stage logs under your configured `work_dir`:

```bash
find wannier_run -maxdepth 2 -type f \
  \( -name '*.out' -o -name '*.wout' \)
```

| File | Information |
| --- | --- |
| `bands_NNNN/projwfc.out` | Projection progress; `NNNN` is the band count |
| `bands_NNNN/fatbands.png`, `.csv`, `_summary.json` | Orbital analysis and numerical summaries |
| `attempt_NNN/preprocess.out` | Wannier preprocessing |
| `attempt_NNN/pw2wannier.out` | Interface progress/errors |
| `attempt_NNN/wannier.wout` | Detailed disentanglement and spread convergence for the default seed |
| `attempt_NNN/wannier.out` | Captured localization console output |
| `decision_NNN_request.json`, `_response.json` | Jev requests and decisions |
| `history.json` | Completed attempts, repairs and band-expansion actions |
| `provenance.json` | Input/code-rule provenance, configuration and actual scratch path |
| `status.json` | Final outcome or a handled failure |

For example, follow the appropriate attempt log:

```bash
tail -F wannier_run/attempt_001/wannier.wout
python3 -m json.tool wannier_run/history.json
python3 -m json.tool wannier_run/status.json
```

Use the latest actual attempt directory; numbering can skip when a decision expands bands. A custom `seed` changes the `wannier.*` filenames. `status.json` need not exist while the first stage is running; history updates after attempts finish.

Success is `"status": "localized"`. The status file identifies the successful attempt containing the checkpoint, Hamiltonian (`*_hr.dat`) and final plan. `analysis_complete` means the requested analysis-only run finished; it does not mean localization occurred. Other outcomes include `stopped`, `attempt_limit`, and `error`.

Exit code 0 means localization or requested analysis completed, 2 means Jev stopped or the attempt budget was exhausted, and 1 means a handled input/API/process error. Read the console log for failures occurring before result initialization.

## Timestamps and projection output

Runner progress messages use the execution machine's local time zone, for example
`[2026-11-04T12:30:00+05:30] Starting pw; output: ...`. They are flushed
immediately, including when redirected to a log with `nohup`. Each external
stage records its start, finish, exit code (or error), and elapsed seconds in
`<work_dir>/stage_events.jsonl`. Jev requests, retries, selections, and Wannier
attempt histories also include timestamps. Raw QE/Wannier output files retain
their original format for analysis tools. Existing runs need to finish before
updating; timestamps apply when the updated runner is launched.

The fat-band stage disables both `lsym` and `kresolveddos` in `projwfc.in`.
Python reads per-band orbital weights from the projection output, so a separate
broadened PDOS calculation is unnecessary. This also avoids the downstream
`partialdos` phase, where an affected QE/MPI installation can hang after printing
the spilling parameter. A busy MPI rank alone does not demonstrate progress.

Timestamps include the machine's UTC offset (`+05:30` means India time). A `TZ`
environment setting is respected. Old log lines retain their original timestamps.
