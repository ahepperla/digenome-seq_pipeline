# Unified cleavage-calling algorithm

## Scope

`bin/call_cleavage.py` calls strand-aware 5-prime endpoint pileups in two
modes:

- `digenome`: paired forward/reverse endpoints representing DSBs
- `ndigenome`: isolated strand endpoints representing SSBs

The same engine supplies alignment filters, local artifact metrics, matched
controls, VCF annotation, and output schemas for both modes.

## Alignment filters

Primary endpoint counts and depth exclude reads that are unmapped, secondary,
supplementary, QC failed, duplicate marked, or below the configured MAPQ.
Secondary/supplementary endpoint support is counted separately.

With `--keep_multimappers`:

- bwa-mem2 emits secondary alignments with `-a`
- the Digenome and nDigenome minimum MAPQ values become `0`
- the supporting-read mean-MAPQ filter becomes `0`
- fastp low-complexity filtering is disabled
- each read is still counted only at BWA's selected primary placement

Secondary and supplementary alignments remain diagnostic. Counting every
reported placement would inflate support in repetitive regions, but
primary-only counting can divide support among equivalent reference copies.

## Coordinates

The internal endpoint coordinate is the 0-based aligned 5-prime reference
base:

- forward: `reference_start`
- reverse: `reference_end - 1`

Soft and hard clips do not extend the aligned endpoint. Insertions do not
advance the reference coordinate; deletions do.

## Streaming endpoint generation

The coordinate-sorted BAM is processed contig by contig. Forward endpoints are
reduced as ordered runs. Reverse endpoints are held in a bounded heap until
the coordinate sort position advances beyond them. Only endpoint runs meeting
the mode's scan threshold are retained for targeted local analysis.

This avoids a genome-sized endpoint dictionary.

## Parallel chunks

The Nextflow workflow divides mapped BAM sequence into a configurable number
of coordinate-based chunks. BAM index mapped-record totals define each
contig's estimated read density. Without a blacklist, large contigs are split
proportionally by genomic coordinate. With a blacklist, boundaries instead
follow cumulative callable mapped work: masked spans contribute zero work and
cuts advance across callable segments. Small contigs can share a chunk. Each
chunk is processed by an independent one-CPU caller.

Every interval has a nonoverlapping 0-based half-open ownership range and a
larger padded scan range. Padding is derived from the configured artifact,
nDigenome opposite-strand, and Digenome pairing neighborhoods. Chunk callers
scan the padded range, but nDigenome calls are owned by their focal endpoint
and Digenome calls are owned by their forward endpoint. This preserves nearby
strand evidence across boundaries while emitting each call once.

When `--genome_blacklist` is supplied, merged BED intervals are subtracted
from each padded scan range before indexed BAM iteration. Focal nDigenome
endpoints, Digenome forward/reverse endpoints, and nDigenome opposite-strand
classification evidence inside the blacklist are excluded. Chunk plans report
owned, callable, and excluded bases plus estimated callable mapped records.
Ownership still covers every base exactly once, including masked spans, so the
finalizer's gap and overlap checks remain unchanged.

Chunk callers write intermediate JSONL rows without applying control q-values
or final filters. A final SQLite-backed merge orders all rows, calculates
Benjamini-Hochberg q-values across the complete sample, applies filters, and
writes the normal TSV, BED, QC, and MultiQC outputs. Before merging, it
requires interval ownership to be gap-free and nonoverlapping across the full
length of every mapped BAM contig. Candidate totals, final filtering, and
Benjamini-Hochberg q-values therefore remain sample-wide.

## nDigenome mode

Each count-qualified endpoint receives local same-strand depth and endpoint
fraction measurements. The caller measures every opposite-strand endpoint in
the configured window and ranks them by:

1. primary count-and-fraction threshold pass
2. ambiguity threshold pass
3. endpoint count
4. endpoint fraction
5. proximity to the focal endpoint
6. coordinate for deterministic ties

Signal classes:

- `SSB`: no meaningful opposite-strand signal
- `POSSIBLE_DSB`: opposite strand passes the primary count and fraction
- `AMBIGUOUS`: opposite strand passes a weaker count or fraction threshold

Only unfiltered `SSB` rows are high confidence.

Default nDigenome criteria:

| Setting | Default | Comparison |
| --- | ---: | --- |
| Focal endpoint count | 10 | `>=` |
| Focal endpoint fraction | 0.20 | `>=` |
| Minimum MAPQ | 1 | `>=` |
| Opposite-strand window | 5 bp | plus or minus |
| Ambiguous opposite count | 3 | `>=` |
| Ambiguous opposite fraction | 0.05 | `>=` |

Classification uses:

```text
POSSIBLE_DSB:
    opposite_count >= 10
    AND opposite_fraction >= 0.20

AMBIGUOUS:
    opposite_count >= 3
    OR opposite_fraction >= 0.05

SSB:
    neither opposite-strand condition is met
```

The Nextflow parameters are `ndigenome_min_count`,
`ndigenome_min_fraction`, `ndigenome_min_mapq`,
`ndigenome_opposite_window`, `ndigenome_ambiguous_min_count`, and
`ndigenome_ambiguous_min_fraction`.

## Digenome mode

Forward and reverse endpoints are eligible to pair when:

```text
abs(reverse_position - (forward_position - overhang)) <= pair_window
```

Candidate pairs are resolved with deterministic one-to-one bipartite
matching. The matching maximizes, in order:

1. pairs that pass every caller threshold
2. total selected pairs
3. total `digenome_pair_score`
4. stable coordinate order

This prevents a high-scoring local pair from consuming endpoints that could
form a larger set of nonconflicting calls. Disconnected endpoint neighborhoods
are solved independently, keeping matching work local to the configured pair
window.

Using endpoint fractions on a zero-to-one scale:

```text
digenome_pair_score = forward_fraction * reverse_fraction
                       * (forward_count + reverse_count) / 4
```

The count, depth, fraction, and pair-score checks use strict greater-than
semantics. Evaluated pairs that fail a threshold remain in the audit TSV.
Only unfiltered `DSB` rows are high confidence.

For historical comparison, `rgen_digenome_score` reproduces the standalone
CRISPR RGEN Tools v1.0 calculation. For each selected forward endpoint it:

1. anchors the historical reverse coordinate at
   `forward_position + overhang - 1`
2. evaluates offsets `-2` through `+2` around both anchors
3. uses unstranded total depth at each coordinate
4. subtracts one from each endpoint count
5. adds the forward-anchored and reverse-anchored contributions using
   single-precision arithmetic

Each directional contribution is:

```text
((count_1 - 1) / total_depth_1)
* ((count_2 - 1) / total_depth_2)
* (count_1 + count_2 - 2)
```

This comparison score follows the historical executable's alignment filter,
including supplementary but excluding unmapped, secondary, QC-failed,
duplicate, and low-MAPQ alignments. It does not affect pairing or filtering.
Blacklist-masked neighborhood positions contribute nothing.

Default Digenome criteria:

| Setting | Default | Comparison |
| --- | ---: | --- |
| Overhang | 0 bp | coordinate adjustment |
| Pair window | 2 bp | plus or minus |
| Minimum MAPQ | 1 | `>=` |
| Forward endpoint count cutoff | 5 | `>` |
| Reverse endpoint count cutoff | 5 | `>` |
| Per-strand depth cutoff | 10 | `>` |
| Per-strand endpoint fraction cutoff | 0.20 | `>` |
| Digenome pair-score cutoff | 2.5 | `>` |

For example, a forward cutoff of `5` requires at least six forward endpoints.
The Nextflow parameters use the `digenome_*` prefix.

## Local artifact metrics

Targeted indexed-BAM fetches calculate:

- same-strand depth and endpoint fraction
- supporting and local MAPQ
- supporting and local NM burden
- 5-prime clipping
- nearby insertion/deletion support
- secondary/supplementary endpoint support
- optional indexed-VCF indel overlap

Nearby indel and known-variant evidence is explicitly reported because
alignment differences from the reference can create artificial endpoint
pileups.

The VCF header must contain every analyzed BAM contig. Declared contig lengths
must match. This turns reference-build and `chr1`/`1` naming mismatches into
clear errors rather than silently missing annotations.

## Matched controls

For each treated candidate, the same coordinate or coordinate pair is measured
in the matched control. The caller reports control counts, depths, fractions,
pseudocount-stabilized fold enrichment, a two-sided Fisher exact p-value, and
a Benjamini-Hochberg q-value.

Controls are optional. Rows without a control are labeled `UNCONTROLLED`.
When a control is present but local depth is below `control_min_depth`, the
row is labeled `INSUFFICIENT_CONTROL_COVERAGE`; fold, p, and q values remain
blank and the row is filtered.

For adequate control depth, fold enrichment uses a 0.5 endpoint pseudocount:

```text
treated_rate = (treated_endpoint_count + 0.5) / (treated_depth + 1)
control_rate = (control_endpoint_count + 0.5) / (control_depth + 1)
fold_enrichment = treated_rate / control_rate
```

Fisher exact p-values compare endpoint and non-endpoint reads. The finalizer
applies Benjamini-Hochberg correction across all candidates for one sample,
after all chunks have been merged.

## Filtering and auditability

Shared artifact reasons include:

- `HIGH_5P_SOFTCLIP`
- `NEARBY_INDEL`
- `KNOWN_INDEL`
- `LOW_SUPPORT_MAPQ`
- `INSUFFICIENT_CONTROL_COVERAGE`
- `HIGH_CONTROL_FRACTION`
- `LOW_CONTROL_FOLD`
- `CONTROL_Q_FAIL`

Digenome threshold reasons are also retained in `filter_reasons`.
nDigenome `POSSIBLE_DSB` and `AMBIGUOUS` rows are filtered from the SSB
high-confidence set.

`signal_classification` records the strand interpretation.
`classification` becomes `ARTIFACT_RISK` when shared artifact evidence is
present. Every evaluated row remains in `*.all.tsv`. Passing rows are also
written to `*.high_confidence.tsv`; filtered `ARTIFACT_RISK` rows go to
`*.artifact.tsv`; all other filtered rows go to `*.manual_review.tsv` for
user review.

Default shared filters:

| Filter | Default failure condition |
| --- | --- |
| Supporting-read mean MAPQ | `< 10` |
| 5-prime soft-clipped fraction | `>= 0.20` |
| Nearby-indel read fraction | `>= 0.20` |
| Known-indel VCF overlap | present |
| Matched-control local depth | `< 1` |
| Matched-control endpoint fraction | `> 0.05` |
| Matched-control fold enrichment | `< 5.0` |
| Matched-control q-value | `> 0.05` |

The three control enrichment conditions are evaluated independently, so a row
can report one, two, or all three detailed reasons. The local artifact window
is 10 bp by default. These settings use the `cleavage_*` Nextflow parameter
prefix.

## Cutoff provenance

The defaults have three different origins:

1. **Digenome reference settings.** The count, depth, fraction, and numeric
   score defaults originate from:

   ```text
   digenome -G 0 -q 1 -f 5 -r 5 -d 10 -R 0.2 -s 2.5
   ```

   The command behavior came from the
   [original Digenome distribution](http://www.rgenome.net/static/digenome-js/digenome).
   The pipeline applies `2.5` to `digenome_pair_score`, so it is not
   numerically equivalent to the standalone command's `-s 2.5` filter.
   The separately reported `rgen_digenome_score` was cross-checked directly
   against the standalone v1.0 executable and is comparison-only.

2. **Published nDigenome focal thresholds.** The defaults of at least 10 reads
   sharing a 5-prime endpoint and at least 20% local fraction come from Kim et
   al., *Unbiased investigation of specificities of prime editing systems in
   human cells*, Nucleic Acids Research (2020),
   [doi:10.1093/nar/gkaa764](https://doi.org/10.1093/nar/gkaa764),
   [PMCID: PMC7544197](https://pmc.ncbi.nlm.nih.gov/articles/PMC7544197/).

3. **Pipeline-selected safeguards.** The MAPQ floor, opposite-strand window,
   ambiguity thresholds, artifact window, clipping and indel limits, and
   control filters are implementation defaults. They are not claimed to come
   from the nDigenome publication.

Every threshold is configurable. Production settings should be calibrated
against known on-target sites, validated off-target sites, matched negative
controls, sequencing depth, library preparation, and the expected nuclease.
