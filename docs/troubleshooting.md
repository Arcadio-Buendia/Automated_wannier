# Troubleshooting and maintenance

[Documentation index](README.md) · [Project README](../README.md)

## Common problems, stopping and cleanup
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

