# Longleaf benchmark: old versus new caller

Run this on a production BAM to check the rebuilt caller against the
pre-refactor caller (commit `b9455e2`) and to confirm the pair-score cutoff.
It runs both callers inside the cleavage image on the same BAM, with the same
chunk count and the pipeline's default thresholds, one step at a time. The old
caller gets the pipeline's `digenome_pair_score_cutoff` (1.1) in place of its
own default of 2.5.

It prints wall time, CPU time, and peak memory for every step, then checks
that every output value matches once the old columns are mapped to the new
layout (`legacy_rows_in_new_layout` in `benchmark.py`). It exits 1 if any
output differs.

```bash
cd /path/to/digenome-seq_pipeline
mkdir -p bench/legacy_src
git archive b9455e2 bin | tar -x -C bench/legacy_src

sbatch -n 1 --mem=32g -t 48:00:00 -o bench/digenome.log --wrap "\
  apptainer exec --bind /proj --bind /work --bind /users \
    containers/cleavage_pysam_v0.23.3.sif \
    python3 tests/benchmark/benchmark.py \
      --legacy-bin bench/legacy_src/bin \
      --analysis digenome \
      --bam /path/to/Sample.sorted.markdup.bam \
      --control-bam /path/to/Control.sorted.markdup.bam \
      --chunks 8 \
      --out bench/digenome"
```

Repeat with `--analysis ndigenome` and `--out bench/ndigenome`. Add `--vcf`,
`--blacklist`, or `--keep-multimappers` to match how the sample is normally
run. Timings are also written to `<out>/timings.json`.

Because the steps run one after another, compare step by step: the slowest
chunk bounds a real run's calling time, and the finalize step now also does
the Digenome pairing.

## Pair-score cutoff

For Digenome, the log also has a table of pair score > 1.1 against
`rgen_digenome_score` > 2.5, the standalone RGEN cutoff, over the pairs that
pass the count, depth, and fraction cutoffs. It gives the agreement at 1.1
and the range of cutoffs from 0.5 to 3.0 with the best agreement. If that
range is far from 1.1, bring the table to the project lead: 1.1 came from a
model of cut sites (`docs/decisions.md`, "The Digenome pair-score cutoff
defaults to 1.1"), and this is its check on real pairs.
