# Calculation configuration

[Documentation index](README.md) · [Project README](../README.md)

## Reuse an existing NSCF calculation
Copy and edit the minimal reuse example:

```bash
cp ~/tools/Automated_wannier/examples/wannier_config.reuse.example.json wannier_config.json
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

## Generate NSCF from an existing SCF calculation
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

## Set spin, target energies, paths and scratch
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

## Multiple cores
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

## Serial Wannier90 builds

A serial Wannier90 executable must run as one process: do not wrap it in a
multiple-rank MPI command. Its output identifies a serial build; multiple copies
can collide on `.werr` and other files. MPI-enabled QE steps can still use multiple
ranks independently. Workflow failures and stops are also reported to the console.
