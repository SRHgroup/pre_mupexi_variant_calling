# Pipeline Rules

These rules are project requirements, not preferences. Follow them for every new pipeline step and every refactor.

## Always Submit Work With `qsub`

- Every real pipeline step must submit a PBS `qsub` job.
- Do not run patient/cohort data processing directly on the login node from `run_pipeline.sh`, `make`, or a step wrapper.
- If a step currently runs work inline, treat that as a bug to fix before expanding the step.
- The wrapper may do lightweight argument parsing, config loading, input prechecks, output-exists checks, and qsub duplicate-submission checks locally.
- Heavy work starts only inside the generated qsub runscript.

## Allowed Local Commands

- `help`, `show-config`, `check`, `check-step`, `watch-step`, and dry-run/precheck commands may run locally.
- Tiny metadata inspection is allowed locally only when it avoids submitting a broken job.
- Debug-only local execution is allowed only when explicitly requested for debugging.

## Required Wrapper Pattern

- `examples/run_pipeline.sh` is the user entrypoint and must call a Make target.
- The Make target must call a step wrapper script.
- The step wrapper must generate a per-step/per-patient runscript under a logs directory.
- The step wrapper must submit that runscript with `qsub`.
- The wrapper should print the submitted job id and where logs/reports will be written.
- Use one job per patient/sample by default.

## Modules And Environment

- Load software modules inside the generated qsub runscript, not only in the interactive shell.
- Step-specific module variables belong in `pipeline_defaults/toolchain.defaults.sh` and can be overridden in `CONFIG`.
- For splicing Python steps, use:

```bash
splicing_module_prereq="tools"
splicing_python_modules="anaconda3/2025.06-1"
```

## Paths And Config

- Do not overload path variables for different concepts.
- BAM-adjacent STAR outputs and splice-junction report outputs may live in different places.
- Splicing STAR junction reports should use `splicing_sjdir`.
- Splicing-derived outputs should use `splicing_outdir` or default to `${datadir}/splicing`.
- Use standardized cohort tags such as `normal_tag` and `tumor_tag`; do not introduce step-specific tumor spellings.
- DNA/RNA labels should be composed from the same tag, e.g. `DNA_${tumor_tag}` and `RNA_${tumor_tag}`.
- Patient-level spl2/spl3/spl3.5/spl4 files should live under `${splicing_outdir}/${patient_id}_RNA_${tumor_tag}`, defaulting `tumor_tag` to `TUMOR`.
- Optional normal-junction filtering should run as `spl3.5`, consuming a compact normal reference from `splicing_normal_junction_ref` or `--normal-ref` and writing a filtered cancer-unique TSV before spl4.
- Normal-reference build and liftover jobs should use the splicing reference qsub wrappers rather than running multi-GB conversions on the login node.
- If both generic and patient-prefixed STAR/splicing files exist for the same sample, keep the patient-prefixed source and skip the generic duplicate.
- Keep cohort-specific paths in `CONFIG`, not hard-coded in scripts.

## Output Behavior

- Steps should skip existing outputs unless `-f` or `FORCE=1` is set.
- Steps should validate required inputs before submission.
- Outputs and logs should be deterministic and discoverable from `CONFIG`, patient id, and step name.

## Current Splicing Reminder

- `spl1`, `spl2`, `spl3`, optional `spl3.5`, and `spl4` must follow the qsub pattern.
- `spl1` should submit the STAR junction shard merge as a qsub job.
- `spl2` should submit the novel junction calling/GTF parsing as a qsub job.
- `spl3` should submit the SSNIP-style event classification/GTF parsing as a qsub job.
- `spl3.5` should submit normal-junction filtering as a qsub job.
- `spl4` should submit neojunction nucleotide/protein sequence reconstruction as a qsub job.
