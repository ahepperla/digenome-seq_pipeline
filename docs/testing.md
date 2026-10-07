# Testing

## Unit and integration-style tests

Run:

```bash
python3 -m pip install -r requirements-test.txt
./tests/run_tests.sh
```

The suite uses Python `unittest` and synthetic indexed BAMs generated with
pysam. It does not require production FASTQs or references. When Nextflow is
installed, it also runs `nextflow lint` on `main.nf` and stub runs of both
calling modes.

Coverage includes:

- samplesheet validation (columns, duplicates, R1/R2, paired-end enforcement)
- sample and control names (character restrictions, control row validation)
- metadata and control consistency
- forward/reverse and complex-CIGAR endpoint coordinates
- duplicate, secondary, supplementary, and MAPQ filtering
- Digenome DSB pairing, overhangs, scores, and one-strand rejection
- deterministic maximum-cardinality Digenome pairing
- nDigenome SSB, possible DSB, and ambiguous calls
- opposite-strand ranking by threshold pass and fraction
- shared clipping, indel, VCF, control, and artifact-risk behavior
- detailed matched-control filter reasons and candidate output tiers
- VCF contig compatibility and insufficient-control-coverage handling
- complete, gap-free, nonoverlapping chunk interval validation
- blacklist-adjusted callable-work balancing, including fully masked contigs
- optional BED/BED.gz blacklist parsing, validation, provenance, and scanning
- Fisher exact and Benjamini-Hochberg calculations
- parameter schema compliance
- complete index reuse
- FASTA fingerprint changes
- partial index quarantine
- stale local lock recovery
- shared cache and lock permissions under a restrictive user umask
- static workflow and container contracts

## Smoke test

Generate fixture inputs:

```bash
python3 tests/fixtures/build_smoke_fixture.py
```

Build the unified image, then run on a local workstation with Apptainer and
Nextflow:

```bash
./containers/build_cleavage.sh

nextflow run . \
  -profile apptainer \
  -c tests/fixtures/smoke.config \
  --input tests/fixtures/tiny_samplesheet.csv \
  --genome tiny \
  --analysis ndigenome \
  --outdir smoke_results \
  -work-dir smoke_work
```

Repeat with:

```bash
nextflow run . \
  -profile apptainer \
  -c tests/fixtures/smoke.config \
  --input tests/fixtures/tiny_samplesheet.csv \
  --genome tiny \
  --analysis digenome \
  --outdir smoke_results_dsb \
  -work-dir smoke_work_dsb
```

The fixture has a DSB at 250 and an SSB at 300, each with 11 forward and 11 reverse
endpoints (or 11 forward-only for the SSB). After Digenome smoke completion,
`Tiny.digenome.all.tsv` should contain at least one data row for the paired
endpoints. After nDigenome smoke completion, the audit output should contain
evidence of the DSB (as `POSSIBLE_DSB` or `AMBIGUOUS` rows) and a passing `SSB`
row at 300. The fixture validates orchestration and basic caller integration,
not biological sensitivity.

## Longleaf validation

On Longleaf, use the SLURM-backed `longleaf` profile rather than the local
`apptainer` profile:

```bash
nextflow run . \
  -profile longleaf \
  -c tests/fixtures/smoke.config \
  --input tests/fixtures/tiny_samplesheet.csv \
  --genome tiny \
  --analysis ndigenome \
  --outdir smoke_results \
  -work-dir smoke_work
```

Repeat with `--analysis digenome`, a separate output directory, and a separate
work directory. The smoke configuration reduces resources only for the tiny
fixture; production runs retain `conf/base.config` resources.

The following cannot be proven by local tests:

- SLURM submission and accounting
- Longleaf paths and mounts (`/proj`, `/work`, `/users`, `/overflow`, `/nas`)
- production SIF execution
- full-depth runtime on human whole-genome sequencing
- biological sensitivity and specificity

Run both smoke modes on Longleaf, then validate with known-positive material,
untreated controls, and representative library preparations. Include a
repetitive-region truth set when validating `--keep_multimappers`. See
`tests/benchmark/README.md` for the Longleaf benchmark kit.
