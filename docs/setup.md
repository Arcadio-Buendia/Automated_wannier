# Installation and input data

[Documentation index](README.md) · [Project README](../README.md)

## Prerequisites
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

## Get the code and prepare Python
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

## Prepare the calculation and API key
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

