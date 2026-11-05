# Automated Wannier

A Jev-directed workflow that takes completed Quantum ESPRESSO SCF or full-mesh NSCF data through orbital analysis and Wannier90 localization.

## Quick start

Requires Python 3.10+, QE (`pw.x`, `projwfc.x`, `pw2wannier90.x`), Wannier90, and a TypeSafe API key. Keep the complete QE save data and pseudopotentials accessible. See [setup](docs/setup.md) for installation details.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cp examples/wannier_config.reuse.example.json wannier_config.json
cp .env.example .env
chmod 600 .env
```

Edit `wannier_config.json` to select your completed NSCF data and a fresh `work_dir`. Put your TypeSafe key in `.env`. For SCF input, use [the fresh-NSCF example](examples/wannier_config.example.json) instead. [Configuration guide](docs/configuration.md).

Run from this folder:

```bash
./run_autonomous.sh
```

To keep running after closing the terminal:

```bash
nohup ./run_autonomous.sh > wannier.log 2>&1 < /dev/null &
disown
```

Follow progress with `tail -F wannier.log`. Timestamps use the machine's local time. The final outcome is in `<work_dir>/status.json`. Use one process for a serial Wannier90 build; configure MPI separately for supported executables.

## Documentation

[Full documentation](docs/README.md) covers setup, configuration, running, troubleshooting, and decision rules. The active Jev rules are in [knowledge/wannier_rules.json](knowledge/wannier_rules.json).

Numerical convergence alone does not establish that the model is scientifically adequate.
