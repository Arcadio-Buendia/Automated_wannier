# Automated Wannier

A runnable, Jev-directed Quantum ESPRESSO → fat-band analysis → Wannier90 workflow. The active decision rules are in `wannier_rules.json`.

## Run from an arbitrary calculation folder

Install Python dependencies from `requirements.txt`. Put `pw.x`, `projwfc.x`, `pw2wannier90.x`, and `wannier90.x` on PATH, or specify command arrays in the configuration. Use compatible modern QE/Wannier90 builds: the baseline targets Wannier90 3.1 or newer with a QE interface supporting SCDM. Pseudo-atomic projection support varies with QE version; an unsupported method triggers a fallback.

The starting calculation needs a converged QE output, its original QE input, the corresponding `prefix.save/data-file-schema.xml` and saved charge density, plus accessible pseudopotentials. A text output alone cannot supply wavefunctions or orbital-resolved bands. The runner creates the required NSCF wavefunctions itself.

Copy `wannier_config.example.json` into your calculation folder as `wannier_config.json` and adjust `qe_input`, `dft_output`, and `save_dir` to your files. `save_dir` can be omitted to derive it from the input's outdir/prefix. All relative paths resolve against the configuration folder. Put `TYPESAFE_API_KEY` in the environment or `.env` beside that configuration.

```bash
python3 -m pip install -r "/path/to/Automated_wannier/requirements.txt"
bash "/path/to/Automated_wannier/run_autonomous.sh" --config wannier_config.json
```

To generate the NSCF/projection results and Python fat-band plots without querying Jev:

```bash
bash "/path/to/Automated_wannier/run_autonomous.sh" --config wannier_config.json --analyze-only
```

Background execution uses `launch_autonomous.sh` with the same arguments. `QE_ENV_FILE` can name an environment setup script. For MPI, use argument arrays such as `"pw": ["mpirun", "-np", "8", "pw.x"]`. Set separate command arrays for the postprocessors as appropriate to your installation. The Wannier command must support standalone serial preprocessing and localization; use its supported launcher if supplying a parallel build.

For collinear magnetic data, select `spin_component: "up"` or `"down"` and use separate work directories for the channels. For SOC/noncollinear data, use `"none"`; spinor multiplicity is read from the saved calculation. `mesh` controls the full unshifted mesh; choose dimensions suitable for the actual geometry, e.g. `[1,1,12]` for a one-dimensional periodic direction. The example `[4,4,4]` is a configurable starting point, not a material-specific recommendation.

## Reuse a completed NSCF calculation

Use `wannier_config.reuse.example.json`. Set `existing_nscf.folder` to the directory where the completed NSCF was run. The defaults expect `nscf.in` and `nscf.out` there; change `input` and `output` if the filenames differ. The save directory is derived from that NSCF input's outdir/prefix, relative to the specified folder. Set `existing_nscf.save_dir` explicitly if the files were relocated; it must point to the actual `prefix.save` directory.

```json
{
  "existing_nscf": {
    "folder": "/path/to/completed_nscf",
    "input": "nscf.in",
    "output": "nscf.out"
  },
  "work_dir": "wannier_reuse_run"
}
```

The reuse example selects the source NSCF folder, its filenames, a new results directory, and a scratch location. It does not need `qe_input`, `dft_output`, `initial_nbnd`, `mesh`, or executable command overrides. Set the placeholder folder and filenames to your actual files. Add `commands` only to override executable paths or enable MPI; the separate parallel example demonstrates those optional settings. The omitted seed, spin selection and target range use the documented defaults; collinear magnetic data still need an explicit spin channel.

The normal launcher command is unchanged. This mode skips the initial pw.x run and starts at projwfc/Python analysis. The full saved NSCF wavefunctions (`wfc*.dat`/HDF5) must be inside the `.save` directory; the program copies it into private scratch so the source stays untouched. External/distributed wavefunction layouts are not automatically collected.

Omit `mesh` to infer the full grid from the saved k points, or specify its matching dimensions. Uniform shifted meshes and arbitrary saved k-point order are retained. A symmetry-reduced grid or a high-symmetry band path is not a substitute for the full Wannier mesh. The automatic grid in the NSCF input is also checked when available. Omit `initial_nbnd`: the saved band count is used. Jev may request more bands later, in which case only that expansion runs a new NSCF in private scratch. Pseudopotentials and pw.x must be available for expansion. Existing projection outputs are not reused; projwfc is rerun on the private copy.

In reuse mode, the NSCF input/output fields above replace the top-level `qe_input`, `dft_output` and `save_dir`. The target reference is the saved NSCF Fermi energy unless you set `fermi_ev` explicitly. Choose a fresh `work_dir` for each trial.

## What the workflow does

1. Copies the original QE save directory into private scratch, retaining the source data. It uses the saved final cell/positions, including for a relaxed input.
2. Generates a full uniform-mesh NSCF input, computes additional bands, and runs `projwfc.x`.
3. Parses orbital weights in Python and writes `fatbands.csv`, `fatbands.png`, and `fatbands_summary.json`. The plot uses uniform-mesh k-point indices, not a high-symmetry path. QE stdout projection weights are rounded; this is sufficient for heuristic shell ranking but is not an exact projectability dataset.
4. Offers Jev concrete, validated plans choosing initialization, shell set, `num_wann`, retained `num_bands`, excluded bands, and energy windows. Jev can request a larger NSCF band pool. Initial orbital choices cover different fractions of the observed target weight. Manual choices preserve complete angular-momentum shells; crystal-field subsets, bond-centered hybrids and custom external projectors are not inferred automatically.
5. Generates `.win` and `.pw2wan`, runs Wannier90 preprocessing, computes interface matrices, and runs localization. It checks generated matrix dimensions before localization.
6. On failure, offers applicable window, orbital, SCDM, or minimizer repairs before method switches. Repairs have finite per-method budgets; failed plans are not repeated. Compatible matrices are reused for window/minimizer changes. Band-pool or initialization changes rebuild the dependent matrices. Elongated-cell neighbor-shell failures can increase the search before giving up; B1 consistency checks remain enabled.
7. Stops successfully only on explicit localization convergence, disentanglement convergence where required, finite final spreads, and nonempty checkpoint/Hamiltonian files. A training-mesh band-energy comparison is reported without imposing a strict accuracy threshold.

The default target is `[-2,2]` eV relative to the source Fermi level. Set `target_relative_ev` for the intended manifold; optionally set `fermi_ev` explicitly. Input windows and SCDM parameters are always written in absolute eV. The source Fermi reference remains fixed across the generated NSCF stages.

## Rules and method selection

- **Bloch phases:** offered only for a sampled isolated band slice with `num_bands == num_wann`; excluded bands remove the rest. Never used with disentanglement. See the [Wannier forum explanation](https://lists.quantum-espresso.org/pipermail/wannier/2021-January/001846.html).
- **SCDM-k:** erfc selection targets occupied/near-Fermi manifolds; Gaussian selection is another choice for bounded energy manifolds. Python attempts a projectability-based erfc fit, reports its quality, and uses a target-aware initial cutoff. Jev can tune the selection. Localization still runs after SCDM initialization. See [tutorial 27](https://wannier90.readthedocs.io/en/stable/tutorials/tutorial_27/) and [tutorial 32](https://wannier90.readthedocs.io/en/stable/tutorials/tutorial_32/).
- **Manual projections:** chosen from fat-band atom/angular-momentum character, with per-site centers and the appropriate spin components. Repairs can change shell coverage or radial extent. See the [projection guide](https://wannier90.readthedocs.io/en/stable/user_guide/wannier90/projections/).
- **Pseudo-atomic projections:** use the orbitals carried by the pseudopotentials, with explicit exclusions and matching projector/Wannier counts. See [tutorial 35](https://wannier90.readthedocs.io/en/stable/tutorials/tutorial_35/).
- **Window repair:** enforce enough outer-window states and at most `num_wann` frozen states at every k-point. Frozen constraints can shrink or be released; the original target remains in the reports. See the [frozen-window forum discussion](https://lists.quantum-espresso.org/pipermail/wannier/2017-February/001303.html).
- **Optional projectability disentanglement:** set `projectability_disentanglement: true` only for a build implementing the relevant keywords. This uses the pseudo-atomic AMN, not rounded plotting weights. It is disabled by default for compatibility.

All source links and decision/action definitions are included in `wannier_rules.json`. The finite candidate ordering and target-energy heuristics are workflow choices; no single method or window is guaranteed to work for every material.

## Outputs, recovery, and validation

The configured work directory contains band-analysis folders, one folder per localization attempt, Jev request/response snapshots, `history.json`, `provenance.json`, and `status.json`. Snapshots contain calculation data sent to TypeSafe; API keys are never written into them. `status.json` names the successful attempt folder and its plan.

Existing runs are preserved: choose a fresh `work_dir` to run again. A lock prevents concurrent writes. After an unclean kill, verify the runner has stopped before removing `.runner.lock`. The original save directory must remain stable while it is copied. Private save copies can consume substantial disk space.

Exit code 0 means localization (or requested analysis-only mode) completed; code 2 means the attempt budget was exhausted or Jev stopped; code 1 means a prerequisite, process, API, or input error. API outages stop the workflow instead of fabricating a Jev decision. JSON limits cap attempts, band counts, and iterations.

Run the contract tests with:

```bash
python3 -m unittest discover -s "/path/to/Automated_wannier/tests" -v
```

Tests use simulated QE/Wannier programs and a local mock Jev endpoint. No real QE/Wannier binaries or material inputs are present in this workspace, so no real-material localization is claimed. Spread convergence is the practical stopping goal; independent-path band accuracy, k-mesh convergence and physical/topological suitability need separate assessment.

## Changes made

Replaced the copied DFT-only entry point with a Wannier pipeline; added Python fat-band analysis, source-linked Wannier JSON rules, a bounded Jev action loop, input generation for each method, matrix invalidation/reuse, spin handling, run provenance and executable integration tests. The original Self Healing DFT repository is untouched.

## Multiple cores

The runner accepts a separate MPI command array for each executable; see `wannier_config.parallel.example.json`. Its ranks (16 for pw.x/pw2wannier, 8 for projwfc, 4 for Wannier90) are starting values to benchmark, not an automatic optimum. Use MPI-enabled executables and a launcher matching the installation or scheduler. Stages run sequentially: the example needs at most 16 MPI ranks at once. Launch the Python/Bash runner once; the command arrays launch MPI inside each stage.

For a pure-MPI starting point, set `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` in the environment before launching. This avoids each rank independently starting a full set of library threads. Hybrid threading can be tuned separately.

NSCF can benefit substantially from MPI; projections, overlap construction and Wannier localization also have parallel implementations. Actual speedup depends on bands, k points, FFT grid, memory traffic and disk I/O. Python fat-band plotting and API requests do not use those MPI ranks. Wannier preprocessing is usually small, and currently uses the same configured Wannier command as localization.

In the [QE 7.5 pw2wannier90 source](https://github.com/QEF/q-e/blob/qe-7.5/PP/src/pw2wannier90.f90), standalone mode supports k-point pools but band groups (`-nb > 1`) are explicitly rejected. For that version, an optional trial is `mpirun -np 16 pw2wannier90.x -nk 4`: four k-point pools with four ranks per pool. Prefer pool counts dividing both the MPI rank count and the relevant k-point count. Older QE versions can reject pools; the provided example omits `-nk` for a simpler baseline. Do not copy all pw.x parallel flags into a postprocessor.

For SCDM, look in pw2wannier output for `Running QRCP in serial` versus `Running QRCP in parallel`. QE 7.5 enables the latter through its ScaLAPACK-QRCP build configuration; an MPI executable alone does not guarantee that factorization is distributed. More ranks can therefore fail to accelerate the dominant SCDM step. Benchmark wall time at a few rank counts on the same calculation before using a whole machine. See the [QE forum performance discussion](https://lists.quantum-espresso.org/pipermail/wannier/2020-November/001823.html) and [Wannier90 MPI usage](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/parameters/).

The same QE source requires an explicit Gamma point for SCDM pivot selection. Shifted NSCF grids without Gamma can still be reused with the other supported initialization methods; SCDM candidates are not offered for those grids. The AMN reader accepts QE's SCDM mu/sigma metadata after the three integer dimensions.

## Configuration selection

The launcher reads `wannier_config.json` in the current working directory by default. Use `--config /path/to/wannier_config.json` when the configuration lives elsewhere. To reuse completed NSCF, the configuration needs a top-level `existing_nscf` entry. To generate NSCF, it needs `qe_input` and `dft_output`. A commands-only parallel configuration does not select a calculation. Missing input settings are rejected before creating a work directory, with the selected configuration path and instructions for both modes.

## Scratch location and foreground launch

Run `./run_autonomous.sh --config /path/to/wannier_config.json` to start in the foreground. Without `--config`, it reads `wannier_config.json` from the current directory. `run_wannier_autonomous.sh` remains available as the original entry point; `launch_autonomous.sh` starts the same workflow in the background.

The `scratch_dir` configuration setting chooses a scratch base directory. It defaults to `./scratch`, relative to the configuration folder. For temporary storage set `"scratch_dir": "/tmp"`. Every invocation creates a separate `wannier-*` directory under that base, so independent calculations do not overwrite each other. The exact directory is printed before the copy starts and recorded as `scratch_dir` in the run's `provenance.json`.

`existing_nscf.save_dir` selects the source data; `scratch_dir` selects where the working copy is created. The source save directory is preserved. Scratch remains after success or failure for inspection; remove the printed per-run directory only after its processes have stopped. Temporary filesystem retention and capacity depend on the remote machine.
