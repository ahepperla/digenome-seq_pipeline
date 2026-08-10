# Testing

## Unit and integration-style tests

Run:

```bash
python3 -m pip install -r requirements-test.txt
./tests/run_tests.sh
```

The suite uses Python `unittest` and synthetic indexed BAMs generated with
pysam. It does not require production FASTQs or references.

Coverage includes:

- original and extended samplesheets
- unknown columns, repeated FASTQs, and identical R1/R2 rejection
- paired-end enforcement
- metadata and control validation
- staged-basename collision detection
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
- chromosome-split boundary padding and single-owner call emission
- blacklist-adjusted callable-work balancing, including fully masked contigs
- optional BED/BED.gz blacklist parsing, validation, provenance, and scanning
- Fisher exact and Benjamini-Hochberg calculations
- parameter schema coverage and preflight path/configuration validation
- complete index reuse
- FASTA fingerprint changes
- partial index quarantine
- stale local lock recovery
- static workflow/container contracts

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

The fixture contains unique paired fragments with:

- 11 forward and 11 reverse endpoints at a synthetic DSB
- 11 forward-only endpoints at a synthetic SSB

After Digenome smoke completion, `Tiny.digenome.all.tsv` should contain at
least one data row. After nDigenome smoke completion, the audit output should
contain the paired-strand evidence plus a passing SSB. The fixture validates
orchestration and basic caller integration, not biological sensitivity.

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

The following cannot be proven by local unit tests:

- SLURM submission and accounting
- visibility and writability of `/proj` paths
- site mount behavior
- production SIF execution
- unified caller runtime on full-depth human WGS
- sensitivity and specificity on known-positive biological material

Run both smoke modes, then known-positive and negative controls before
production use. Include a repetitive-region truth set when validating
`--keep_multimappers`.
