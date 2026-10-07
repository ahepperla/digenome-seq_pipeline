# Testing

## Local suite

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -r requirements-test.txt
export PATH="$PWD/.venv/bin:$PATH"
./tests/run_tests.sh
```

The suite needs Python with pysam 0.23.3 (matching the cleavage image) and no
production data. It runs `bash -n` on the shell scripts and every unittest:

- **Golden outputs** (`test_golden.py`): six synthetic scenarios in both modes
  at 1, 2, and 4 chunks must match `tests/golden/expected/` byte for byte.
  When written, those files were checked value by value against the
  pre-refactor caller (commit `b9455e2`). `python3 tests/golden/build_expected.py`
  rewrites them; do that only for an intended output change.
- **Caller units**: endpoints and read filters, site metrics, the RGEN score,
  nDigenome classes and ranking, Digenome pairing and matching, controls,
  Fisher and Benjamini–Hochberg, filters and tiers, the blacklist, the chunk
  plan, and every finalizer check. Threshold tests sit exactly on each
  boundary, so a `>` turned into `>=` fails.
- **Long reads** (`test_long_reads.py`): both aligned ends count, as `+`/`−`
  by side in Digenome and on the read's strand in nDigenome; a one-base
  alignment is one endpoint; Digenome depth is every covering read, and the
  combined depth counts each read once in treated and control; a pair split
  across chunks; the clip at the endpoint's own end; a nick gives two SSB
  rows and a cut gives POSSIBLE_DSB rows.
- **Samplesheet** checks, one test per rule, for both input formats.
- **Index cache** (`test_index_cache.py`) with a fake `bwa-mem2`: reuse,
  content addressing, version isolation, quarantine, stale-lock recovery,
  permissions under a restrictive umask.
- **Workflow** (`test_workflow.py`): parameter names and defaults agree across
  `nextflow.config`, `nextflow_schema.json`, and `parameters.md`; container
  provenance; the smoke fixture. With `nextflow` installed it also runs
  `nextflow lint` and stub runs of both modes and of `--long_reads`
  (`-stub-run`, no tools needed), checking every published file, the
  processes that ran, one chunk task per chunk per called sample, that
  controls are not called, and the `--keep_multimappers` and `--long_reads`
  settings. The Digenome run names its genome by an alias in another case.
  The long-read run has no `--genome`, and its treated and control BAMs share
  a file name.

The workflow tests use `$NEXTFLOW` if set, else `nextflow` on PATH. Run them
on Longleaf's version with `NEXTFLOW=nextflow-25.04.7 ./tests/run_tests.sh`.

Nextflow 25.10 and later parse command-line parameters as strings under their
strict parser, which nf-schema 2.5.1 then rejects. For local runs with a newer
Nextflow, set `NXF_SYNTAX_PARSER=v1` to match Longleaf's 25.04; the stub tests
do.

## Smoke test with real tools

```bash
python3 tests/fixtures/build_smoke_fixture.py      # writes tests/fixtures/generated/
nextflow run . -profile apptainer,test --analysis ndigenome --outdir smoke_ndigenome -work-dir smoke_work_ndigenome
nextflow run . -profile apptainer,test --analysis digenome --outdir smoke_digenome -work-dir smoke_work_digenome
```

On Longleaf use `-profile longleaf,test`. The `test` profile points at the
generated samplesheet and `tests/fixtures/tiny.fa`, uses three chunks, and caps
each task at 2 CPUs, 4 GB, and 1 hour.

The fixture has 33 uniquely aligned read pairs: a DSB at position 250 (11
forward reads start there and 11 reverse reads end there) and an SSB at 300
(11 forward reads start there). Expect a Digenome pair at 250 and, in
nDigenome mode, a passing SSB at 300 plus POSSIBLE_DSB rows at 250. The
fixture checks orchestration and integration, not biological sensitivity.

To smoke-test `--long_reads`, give the BAMs from the Digenome smoke run as
long-read input. The fixture's reads aren't long reads, so only check that
the run finishes and publishes the cleavage tables:

```bash
printf 'sample,bam\nTiny,%s\n' "$PWD/smoke_digenome/bam/Tiny.sorted.markdup.bam" > smoke_long_reads.csv
nextflow run . -profile apptainer,test --long_reads --input smoke_long_reads.csv --outdir smoke_long_reads -work-dir smoke_work_long_reads
```

## What only Longleaf can show

- SLURM submission and accounting
- Longleaf paths and mounts (`/proj`, `/work`, `/users`, `/overflow`, `/nas`)
- the production images
- full-depth runtime
- biological sensitivity and specificity

Before production use, run both smoke modes on Longleaf, run the benchmark in
[tests/benchmark/README.md](../tests/benchmark/README.md) on a real sample,
then validate with known-positive material, untreated controls, and
representative library preparations. Include a repetitive-region truth set
when validating `--keep_multimappers`.
