# Automated Wannier

A Jev-directed Quantum ESPRESSO → orbital-resolved band analysis → Wannier90 workflow. Start from either a completed SCF calculation or a completed full-mesh NSCF calculation. Python prepares valid candidate settings, Jev selects a method or repair, and the runner executes the resulting Wannierization stages.

The active decision rules and their research sources are in [wannier_rules.json](wannier_rules.json). Configuration examples are provided for a [new NSCF calculation](wannier_config.example.json), [reusing NSCF](wannier_config.reuse.example.json), and [MPI execution](wannier_config.parallel.example.json).

## 1. Prerequisites

The commands below assume a Linux system with Bash and Git. Use a supported Python 3 installation; Python 3.10 or newer is a practical starting point. The Python requirements are NumPy, SciPy, and Matplotlib.

Install or load these external programs before running the workflow:

| Program | Used for |
| --- | --- |
| `pw.x` | Generating NSCF wavefunctions, or expanding an existing band pool |
| `projwfc.x` | Atomic-orbital projections and fat-band analysis |
| `pw2wannier90.x` | Preparing Wannier overlap, projection, and eigenvalue files |
| `wannier90.x` | Preprocessing, disentanglement, and localization |
| An MPI launcher, if using multiple processes | Launching the MPI-enabled programs |

The baseline targets Wannier90 3.1 or newer and a QE interface with SCDM support. Pseudo-atomic projections and projectability disentanglement require corresponding build support. Use mutually compatible QE executables and an MPI launcher matching their installation. This repository's Python requirements do not install QE, Wannier90, or MPI; load them using your machine's installation or cluster modules.

An internet connection and a TypeSafe API key are needed for Jev decisions. `--analyze-only` does not call Jev, but still runs the required numerical analysis stages. The workflow sends orbital summaries, candidate plans, failure excerpts and remediation history to TypeSafe when requesting decisions.

## 2. Get the code and prepare Python

Keep the shared code and your calculation folders separate:

```text
~/tools/Automated_wannier/          scripts, rules, examples
~/calculations/example/            configuration, credentials, input data, results
```

For a public repository, clone over HTTPS:

```bash
mkdir -p ~/tools
git clone https://github.com/Arcadio-Buendia/Automated_wannier.git \
  ~/tools/Automated_wannier

python3 -m venv ~/venvs/wannier
source ~/venvs/wannier/bin/activate
python -m pip install -r ~/tools/Automated_wannier/requirements.txt
```

Activate that environment in each new terminal before running the workflow. Check that the scientific programs are available:

```bash
command -v pw.x projwfc.x pw2wannier90.x wannier90.x
```

If they are not on `PATH`, load the site's software environment or set absolute executable paths in the configuration's `commands` arrays. An environment setup script can be sourced by the launcher:

```bash
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export QE_ENV_FILE="/path/to/software-environment.sh"
```

For non-interactive/background execution, ensure that environment is available to the job too. Relative executable paths are resolved inside stage directories; use `PATH` or absolute paths.

## 3. Prepare the calculation and API key

A text output by itself is insufficient. Keep the corresponding save data and pseudopotentials. For example:

```text
~/calculations/example/
├── wannier_config.json
├── .env
├── nscf_run/
│   ├── nscf.in
│   ├── nscf.out
│   └── output/
│       └── sample.save/
│           ├── data-file-schema.xml
│           ├── charge-density.dat       or the HDF5 equivalent
│           ├── wfc*.dat                 or the HDF5 equivalents
│           └── ...                      retain the complete save directory
└── pseudo/
    └── ...UPF
```

Choose one of the two modes below. The supplied input/output filenames must match your actual files, and the completed output must contain `JOB DONE` without `convergence NOT achieved`.

Create `.env` beside your configuration:

```bash
mkdir -p ~/calculations/example
cd ~/calculations/example
cp ~/tools/Automated_wannier/.env.example .env
chmod 600 .env
nano .env
```

Replace its placeholder with your key:

```text
TYPESAFE_API_KEY=your_actual_key
```

The runner first checks the exported `TYPESAFE_API_KEY` environment variable. If it is absent or empty, it reads `.env` beside the selected configuration. An exported old key therefore takes precedence over a newer `.env` value. Credentials and personal calculation files belong outside the shared code checkout.

## 4A. Reuse an existing NSCF calculation

Copy and edit the minimal reuse example:

```bash
cp ~/tools/Automated_wannier/wannier_config.reuse.example.json wannier_config.json
nano wannier_config.json
```

A complete minimal configuration looks like this:

```json
{
  "existing_nscf": {
    "folder": "/absolute/path/to/nscf_run",
    "input": "nscf.in",
    "output": "nscf.out"
  },
  "work_dir": "wannier_run",
  "scratch_dir": "/tmp"
}
```

Set the source folder and filenames. The source `.save` path is normally derived from `outdir` and `prefix` in the NSCF input. If the save directory was moved, add `"save_dir": "output/sample.save"` **inside `existing_nscf`**, or supply its absolute path.

This mode skips the initial `pw.x` calculation and starts with projections. It requires a complete uniform k-point mesh and saved wavefunctions inside the `.save` directory. A high-symmetry band path or a symmetry-reduced mesh is insufficient. Uniform shifted meshes and arbitrary saved k-point order are supported; the runner checks the grid and keeps its ordering. External/distributed wavefunction layouts are not automatically collected.

Omit `initial_nbnd`: the saved band count is used. Omit `mesh` to infer it, or supply matching dimensions. Jev can request more bands later; that expansion will require `pw.x` and accessible pseudopotentials. Existing `projwfc` outputs are not reused.

## 4B. Generate NSCF from an existing SCF calculation

Use `wannier_config.example.json` instead. You need the original QE input, its completed output, the matching save directory and charge density, and accessible pseudopotentials. Existing NSCF wavefunctions are not required in this mode.

```json
{
  "qe_input": "scf.in",
  "dft_output": "scf.out",
  "save_dir": "output/sample.save",
  "work_dir": "wannier_run",
  "scratch_dir": "./scratch",
  "mesh": [4, 4, 4],
  "spin_component": "none"
}
```

The runner uses the saved final cell and positions, copies the save data, and computes a full unshifted NSCF mesh. Choose `mesh` for the actual periodic directions; `[4,4,4]` is only an example. The QE input's pseudopotential paths must resolve on this machine.

## 5. Set spin, target energies, paths and scratch

| Setting | Meaning |
| --- | --- |
| `spin_component: "none"` | Scalar nonmagnetic or noncollinear/SOC data |
| `spin_component: "up"` or `"down"` | One collinear magnetic spin channel; run the channels separately with different `work_dir` values |
| `target_relative_ev` | Target range relative to the reference Fermi level; default `[-2,2]` eV |
| `fermi_ev` | Optional explicit absolute Fermi-energy reference |
| `work_dir` | New results directory for this invocation |
| `scratch_dir` | Base directory for the working copy; default `./scratch` |
| `commands` | Optional executable/MPI command arrays |

`"0"` is not a valid spin selection. SOC is read from the saved calculation; keep `spin_component` as `"none"` for spinor data. Target energies use `E-Ef`; the generated Wannier windows and SCDM parameters are written in absolute eV.

Top-level paths, including `work_dir` and `scratch_dir`, are relative to the **configuration folder**, regardless of where the scripts are installed. `existing_nscf.folder` is also relative to that folder; its `input`, `output`, and `save_dir` paths are relative to the specified NSCF folder. Absolute paths override those bases. `pseudo_dir` is interpreted relative to the NSCF folder in reuse mode and relative to the configuration folder in SCF mode. Use absolute pseudopotential paths if the original input was stored or run elsewhere.

The source save directory and scratch are different locations: `save_dir` selects existing data; `scratch_dir` selects where the working copy is made. For `"scratch_dir": "/tmp"`, each invocation creates a separate `/tmp/wannier-*` directory with owner-only directory permissions. The exact path is printed before copying and recorded in `provenance.json`. The original save is preserved.

Allow disk space for a full save copy plus generated files. Scratch remains after success or failure. `/tmp` may be cleaned by the machine; copy any results you need before that happens. Choose a persistent scratch filesystem for long jobs where appropriate.

## 6. Launch in the foreground or detach

From the calculation folder, with the Python environment active:

```bash
bash ~/tools/Automated_wannier/run_autonomous.sh --config wannier_config.json
```

If the scripts and configuration are in the current folder, `./run_autonomous.sh` is equivalent. Without `--config`, the launcher reads `wannier_config.json` from the current working directory. `run_wannier_autonomous.sh` remains a compatible entry point.

To keep running after closing the console, start a new job with:

```bash
nohup bash ~/tools/Automated_wannier/run_autonomous.sh \
  --config wannier_config.json > wannier.log 2>&1 < /dev/null &
echo $! > wannier.pid
disown
```

Use `disown` from an interactive Bash shell. Do not start another copy of an already running job. `launch_autonomous.sh --config ...` is an alternative background launcher; it appends console output to `wannier_autonomous.log` in the current directory. On a scheduled cluster, put the launch command in the site's batch job and use its allocation and launcher.

An analysis-only trial uses:

```bash
bash ~/tools/Automated_wannier/run_autonomous.sh \
  --config wannier_config.json --analyze-only
```

This still copies save data, runs NSCF if required, runs `projwfc`, and produces Python analysis. It makes no Jev requests and performs no localization. Use a **fresh `work_dir` for a later full run**; the program does not resume an analysis-only directory.

## 7. Multiple cores

Set a separate command array for each program, for example:

```json
{
  "commands": {
    "pw": ["mpirun", "-np", "16", "pw.x"],
    "projwfc": ["mpirun", "-np", "8", "projwfc.x"],
    "pw2wannier": ["mpirun", "-np", "16", "pw2wannier90.x"],
    "wannier": ["mpirun", "-np", "4", "wannier90.x"]
  }
}
```

Add this block to a valid calculation configuration; a commands-only file does not select input data. The parallel example includes an existing-NSCF selection. These process counts are starting values to benchmark. Use MPI-enabled executables; the stages run sequentially, so this example uses at most 16 ranks at once. Launch the workflow once; its command arrays launch the MPI stages. Wannier preprocessing currently uses the same command as localization.

For an initial pure-MPI run:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

NSCF, projections, overlap construction and localization can benefit from MPI. Scaling depends on the problem, libraries and I/O. In the [QE 7.5 interface source](https://github.com/QEF/q-e/blob/qe-7.5/PP/src/pw2wannier90.f90), standalone mode supports k-point pools but rejects band groups (`-nb > 1`). An optional version-appropriate trial is `mpirun -np 16 pw2wannier90.x -nk 4`; older builds may reject pools. Prefer pool counts dividing both ranks and relevant k-point counts.

SCDM's QR factorization may remain serial. Look for `Running QRCP in serial` or `Running QRCP in parallel` in interface output; distributed QRCP needs appropriate build support. MPI alone does not guarantee acceleration of that step. See the [QE forum performance discussion](https://lists.quantum-espresso.org/pipermail/wannier/2020-November/001823.html) and [Wannier90 MPI instructions](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/parameters/).

## 8. Follow logs and inspect results

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

## 9. Common problems, stopping and cleanup

| Symptom | Action |
| --- | --- |
| Missing `qe_input` / invalid input selection | Add `existing_nscf` for reuse, or both `qe_input` and `dft_output` for new NSCF; select the intended file using `--config` |
| Missing `data-file-schema.xml` | Point `save_dir` to the actual relocated `.save` directory |
| No saved wavefunctions | Restore the full NSCF save data or generate NSCF again; text output alone is insufficient |
| Mesh validation failure | Supply a full uniform mesh with matching dimensions, not an irreducible grid or band path |
| Spin selection error | Use `none`, `up`, or `down` as appropriate; separate collinear channels |
| Executable/library not found | Load the site's environment or set compatible absolute commands and library paths |
| TypeSafe key missing or rejected | Check environment-first key selection and `.env` beside the config; API failures stop rather than invent a decision |
| Existing `history.json` in `work_dir` | Choose a fresh results directory; automatic cross-run resume is not implemented |
| `.runner.lock` exists | Check for the active job first; remove a stale lock only after its workflow and children have stopped |
| Scratch full or removed by the system | Choose a suitable `scratch_dir` and start with a fresh results directory |

For a foreground run, use Ctrl+C to interrupt. For background runs, use the saved PID and `ps` to identify the workflow and its MPI child processes before stopping that job; stopping only a wrapper may leave numerical workers alive. Do not delete scratch until the associated processes have stopped. Keep calculation data separate when replacing or updating a code checkout.

Update the shared code between runs:

```bash
git -C ~/tools/Automated_wannier pull --ff-only
```

Reload Python requirements if they change. Existing runs and source data are preserved; rerunning is an explicit new run with a new `work_dir`.

## 10. Decision rules and configurable limits

`wannier_rules.json` contains the decision instructions, supported method descriptions, actions, defaults, bounds and source links. A different complete rules file can be selected using top-level `rules_file`. The executable handlers support the implemented methods/actions; adding an arbitrary new JSON action does not create a new handler.

| Initialization | Selection constraint |
| --- | --- |
| Bloch phases | Sampled isolated retained band manifold, with `num_bands == num_wann`; never combined with disentanglement |
| SCDM erfc | Occupied/near-Fermi manifolds; target-aware selection and optional projectability fit |
| SCDM Gaussian | Alternative bounded-energy manifold selection |
| Manual projections | Per-site angular-momentum shells ranked from projected weight, with spinor multiplicity |
| Pseudo-atomic projections | Orbitals supplied in pseudopotentials, with exclusions and matching projector count |

The QE SCDM implementation used here requires an explicit Gamma point; shifted meshes without Gamma are offered other methods. At every k-point, an outer window must contain at least `num_wann` states and a frozen window at most `num_wann`. Applicable window, orbital, SCDM and minimizer repairs are offered before another initialization method. Compatible matrices are reused; changed band pools/projections rebuild them. Projectability disentanglement is disabled by default and must only be enabled for a supporting build.

Configuration options include `max_attempts` (default 48), `stage_timeout_seconds` (default 86400), `initial_nbnd` for new NSCF, `model` (default `jev-latest`), and `projectability_disentanglement`. The rules also cap generated band pools and iteration budgets. The API timeout is separate from the numerical-stage timeout. API credentials are not stored in decision snapshots; calculation summaries and errors are.

Python produces fat-band plots along uniform-mesh indices rather than a high-symmetry path. QE stdout weights are rounded and are used for heuristic shell ranking, not as exact normalized AMN projectabilities. Full angular-momentum shells are selected; crystal-field subsets, bond-centered hybrids and custom external projectors are not inferred automatically.

## Validation and scientific scope

The stopping goal is explicit localization convergence, disentanglement convergence where required, finite final spreads, and nonempty checkpoint/Hamiltonian files. A training-mesh energy comparison is reported without imposing a strict accuracy threshold. Independent-path interpolation accuracy, k-mesh convergence and physical/topological suitability require separate assessment. Numerical convergence alone does not establish that the model is scientifically adequate.

The repository's contract tests use simulated QE/Wannier programs and a local mock Jev endpoint. They verify orchestration and recovery behavior, not universal real-material convergence:

```bash
python3 -m unittest discover -s ~/tools/Automated_wannier/tests -v
```

Method references include the [projection guide](https://wannier90.readthedocs.io/en/stable/user_guide/wannier90/projections/), [SCDM tutorial](https://wannier90.readthedocs.io/en/stable/tutorials/tutorial_27/), [projectability fit tutorial](https://wannier90.readthedocs.io/en/stable/tutorials/tutorial_32/), [pseudo-atomic tutorial](https://wannier90.readthedocs.io/en/stable/tutorials/tutorial_35/), [Bloch-phase forum discussion](https://lists.quantum-espresso.org/pipermail/wannier/2021-January/001846.html), and [frozen-window forum discussion](https://lists.quantum-espresso.org/pipermail/wannier/2017-February/001303.html). Additional sources are in the rules file.

Runner progress messages use UTC timestamps, for example
`[2026-11-04T12:30:00+00:00] Starting pw; output: ...`. They are flushed
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
