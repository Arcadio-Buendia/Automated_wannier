# Decision rules and validation

[Documentation index](README.md) · [Project README](../README.md)

## Decision rules and configurable limits
`knowledge/wannier_rules.json` contains the decision instructions, supported method descriptions, actions, defaults, bounds and source links. A different complete rules file can be selected using top-level `rules_file`. The executable handlers support the implemented methods/actions; adding an arbitrary new JSON action does not create a new handler.

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

