# Cleavage calling

`bin/cleavage/` calls strand-aware read-end pileups in two modes:

- `digenome`: forward and reverse endpoints paired into double-strand breaks (DSBs);
- `ndigenome`: isolated strand endpoints, i.e. single-strand breaks or nicks (SSBs).

Both modes share the read filters, local artifact metrics, matched controls,
known-indel annotation, and output tiers described here. Short reads
contribute their 5′ ends; long reads (`--long_reads`) contribute both ends
(see [Long reads](#long-reads)). Every threshold is a pipeline parameter; see
[parameters.md](parameters.md).

## Reads that count

Endpoint counts and depths use primary alignments only. Unmapped, secondary,
supplementary, QC-failed, and duplicate-marked reads are excluded, and so are
reads below the mode's minimum MAPQ (`>=` passes). Secondary and supplementary
alignments ending at a site are reported as `secondary_endpoint_count` for
diagnosis only.

With `--keep_multimappers`, bwa-mem2 runs with `-a`, both minimum MAPQs and the
support mean-MAPQ filter become 0, and fastp's low-complexity filter is off.
Each read still counts once, at BWA's selected primary placement. Counting
every reported placement would inflate support in repeats; counting only the
primary placement can split support among equivalent copies.

## Coordinates

The endpoint of a short read is the 0-based aligned 5′ reference base:
`reference_start` for forward reads and `reference_end - 1` for reverse reads.
Soft and hard clips never extend it. Insertions do not advance the reference
coordinate; deletions do. Tables report both 0-based and 1-based columns; the
BED is 0-based half-open.

## Chunks

Each sample is called in `--cleavage_chunks` coordinate chunks, each a
one-CPU task, followed by one finalize task.

**Plan.** Every chunk task computes the same plan from the BAM index and the
blacklist, using exact rational arithmetic so that all tasks agree. Each
mapped contig's mapped records are assumed spread evenly along it, counting
only bases outside the blacklist. The genome's callable work is cut into
`n` equal parts: large contigs are split at those cuts, small ones share a
chunk, and a chunk can be empty when there is too little data to split. The
owned ranges are 0-based half-open, never overlap, and cover every base of
every mapped contig, blacklisted spans included. The finalizer publishes the
plan as `pipeline_info/cleavage_chunks/<sample>.cleavage_chunks.tsv`.

**Calling.** A chunk scans only the ranges it owns, minus the blacklist, so
each endpoint is found by exactly one chunk. Everything else a call needs
(local metrics, the opposite-strand search, control measurements, RGEN
counts, VCF lookups) is read with its own targeted fetch and is not limited
to the chunk. Chunk files are sorted by contig (BAM header order), position,
and strand. They carry no q-values or filters.

**Finalize.** The finalizer checks that all `n` chunks ran with the same
settings, contigs, plan, and blacklist, that the plan owns every mapped contig
exactly once, that every record lies in its own chunk's ranges, and that no
records are missing or repeated. It then streams the chunk files together,
builds the rows (pairing Digenome endpoints, below), computes
Benjamini–Hochberg q-values over every controlled row in the sample, applies
the filters, and writes the outputs.

Because every endpoint is found once and every decision that spans chunks is
made in the finalizer, the number of chunks changes runtime, never results.

## nDigenome

A focal endpoint is called when at least `ndigenome_min_count` reads end on it
and they are at least `ndigenome_min_fraction` of the same-strand depth there
(both `>=`). The caller then measures every opposite-strand endpoint within
`ndigenome_opposite_window` bases, ignoring blacklisted ones, and keeps the
strongest, ranked by:

1. passes the primary count and fraction;
2. passes the ambiguity count or fraction;
3. endpoint count;
4. endpoint fraction;
5. closeness to the focal endpoint;
6. lower coordinate, for deterministic ties.

That endpoint sets the signal class:

```text
POSSIBLE_DSB: opposite_count >= ndigenome_min_count
              AND opposite_fraction >= ndigenome_min_fraction
AMBIGUOUS:    opposite_count >= ndigenome_ambiguous_min_count
              OR opposite_fraction >= ndigenome_ambiguous_min_fraction
SSB:          neither
```

Only unfiltered SSB rows are high confidence; POSSIBLE_DSB and AMBIGUOUS rows
are filtered with their class as the reason.

## Digenome

A forward endpoint `f` and a reverse endpoint `r` can pair when

```text
abs(r - (f - overhang)) <= pair_window
```

Only endpoints with more reads than the smaller of the forward and reverse
count cutoffs are considered. Each chunk records every owned endpoint that has
at least one possible partner, with its treated and control metrics, nearby
known indels, and (for forward endpoints) the RGEN scores. To find partners
owned by a neighboring chunk, it scans `|overhang| + pair_window` bases past
its own ranges.

The finalizer builds every candidate pair on a contig and selects a
deterministic one-to-one matching. Candidates that share an endpoint, directly
or through a chain, form a component, and each component is solved as a
min-cost flow that maximizes, in order:

1. pairs that pass every caller threshold;
2. the number of pairs;
3. the total `digenome_pair_score`, in exact integer units;
4. stable coordinate order.

This keeps one high-scoring pair from consuming endpoints that could form a
larger set of non-conflicting calls. Only selected pairs are reported;
competing pairs the matching did not choose are not written.

```text
digenome_pair_score = forward_fraction * reverse_fraction
                      * (forward_count + reverse_count) / 4
```

Caller thresholds are strict `>`: forward and reverse counts, each strand's
depth, each strand's fraction, and the pair score. A cutoff of 5 therefore
needs at least 6 reads. Pairs that fail one are reported with reasons such as
`LOW_FORWARD_COUNT` and `LOW_DIGENOME_PAIR_SCORE`. The `combined_*` columns
pool both strands (for long reads, each read once); they feed the control test
and the BED score.

### RGEN comparison score

`rgen_digenome_score` reproduces the standalone CRISPR RGEN Tools v1.0 score,
for comparison with historical results only. It never affects pairing,
filtering, or tiers, and its values are not comparable with pair-score
cutoffs. For each forward endpoint it:

1. anchors the historical reverse coordinate at `forward + overhang - 1`;
2. evaluates offsets −2 to +2 around both anchors;
3. uses unstranded total depth at each base;
4. subtracts one from each endpoint count;
5. adds the forward- and reverse-anchored contributions in single precision:

```text
((count_1 - 1) / depth_1) * ((count_2 - 1) / depth_2) * (count_1 + count_2 - 2)
```

It follows the executable's read filter, which keeps supplementary alignments
but excludes unmapped, secondary, QC-failed, duplicate, and low-MAPQ ones.
Blacklisted and off-contig positions contribute nothing.

## Long reads

With `--long_reads`, each sample is one aligned BAM (ONT or PacBio): no
trimming, alignment, or `--genome`. A BAM that isn't coordinate-sorted with an
index beside it is sorted and indexed first.
A long read spans a whole molecule, so both of its aligned ends are molecule
ends, and both count. The read filters and thresholds are the same, apart
from the artifact limits below.

| Data | Endpoints of one read |
| --- | --- |
| short reads | the 5′ end: `reference_start` (+) or `reference_end - 1` (−) |
| long reads, Digenome | `reference_start` as `+` and `reference_end - 1` as `−`, whatever the read's orientation |
| long reads, nDigenome | both ends, each on the read's own strand |

Clips never extend an end, and an alignment whose two ends fall on one base
contributes one endpoint, its left end (`+` in Digenome). Depth at an endpoint
counts the reads covering it that can end there on that strand: every
covering read for long-read Digenome, and reads on that strand otherwise.

**Digenome.** A cut between `p - 1` and `p` leaves molecules whose left end is
`p` and molecules whose right end is `p - 1`, sequenced in either orientation.
Their ends become `+` endpoints at `p` and `−` endpoints at `p - 1`, which pair
exactly as short-read endpoints do. The pair score uses each side's own count
and depth. For the combined columns, the control fraction, the fold, and the
Fisher test, the combined depth counts each read covering either endpoint
once; summing the two sides' depths would count every spanning read twice.

**nDigenome.** A read's strand is the DNA strand of its molecule. In a library
that keeps strands, a nick gives two same-strand stacks one base apart, so two
SSB rows; a double-strand cut adds ends on the other strand within the window,
and those rows become POSSIBLE_DSB. A standard double-stranded ligation
library leaves no read ends at nicks. Long reads need not be paired.

**Artifact limits.** `--long_reads` sets `cleavage_max_softclip_fraction` and
`cleavage_max_indel_fraction` to 1.0, because nanopore reads routinely carry
small indels and clipped ends. The metrics are still reported, and
`HIGH_5P_SOFTCLIP` still applies when every supporting read is clipped at
that end, for example by adapter remnants; check `softclip_fraction` before
trusting the artifact tier.

**RGEN score.** `rgen_digenome_score` counts both ends of long reads too, and
it keeps supplementary alignments as the standalone tool does, so the split
points of chimeric reads enter it. It remains comparison only.

## Local artifact metrics

Around each endpoint (`cleavage_artifact_window` bases each side) the caller
measures supporting and local MAPQ and NM, the fraction of supporting reads
soft- or hard-clipped at the end that forms the endpoint (the 5′ end for short
reads), the fraction of local alignments with a CIGAR insertion or
deletion touching the window (and the most common such indel), secondary
endpoint support, and indel records from the optional VCF. For a Digenome
pair, the two strands' metrics are pooled; the reported indel comes from the
strand with the higher indel fraction, then read count, and from the forward
strand on a tie.

The VCF must be bgzipped and indexed, declare every analyzed contig, and agree
with the BAM on any declared contig length. A mismatch (another build, or
`chr1` versus `1`) stops the run rather than silently losing annotations.

## Matched controls

The same coordinate, or both coordinates of a pair, is measured in the matched
control. A row without a control is `UNCONTROLLED` and is not filtered for
lacking one. When control depth is below `cleavage_control_min_depth`, the row
is `INSUFFICIENT_CONTROL_COVERAGE`: fold, p, and q stay blank and it is
filtered. Otherwise:

```text
treated_rate    = (treated_endpoint_count + 0.5) / (treated_depth + 1)
control_rate    = (control_endpoint_count + 0.5) / (control_depth + 1)
fold_enrichment = treated_rate / control_rate
```

A two-sided Fisher exact test (implemented here, no SciPy) compares endpoint
and non-endpoint reads in treated and control. Benjamini–Hochberg q-values are
computed in the finalizer across every controlled row of the sample.

## Filters and tiers

Shared artifact reasons, and when they apply:

| Reason | Condition |
| --- | --- |
| `HIGH_5P_SOFTCLIP` | soft-clipped fraction `>=` `cleavage_max_softclip_fraction` |
| `NEARBY_INDEL` | indel fraction `>=` `cleavage_max_indel_fraction` |
| `KNOWN_INDEL` | a VCF indel within the artifact window |
| `LOW_SUPPORT_MAPQ` | supporting mean MAPQ `<` `cleavage_min_support_mean_mapq` |
| `INSUFFICIENT_CONTROL_COVERAGE` | control depth `<` `cleavage_control_min_depth` |
| `HIGH_CONTROL_FRACTION` | control endpoint fraction `>` `cleavage_control_max_fraction` |
| `LOW_CONTROL_FOLD` | fold enrichment `<` `cleavage_control_min_fold` |
| `CONTROL_Q_FAIL` | q-value `>` `cleavage_control_max_q` |

The three control conditions are checked independently, so a row can report
one, two, or all three. `filter_reasons` lists Digenome threshold reasons
first, then artifact reasons, then (nDigenome) a non-SSB signal class.

Each row gets a `tier`:

- `high_confidence`: no reasons (`filter_status` PASS); also written to the BED;
- `artifact`: at least one artifact reason;
- `manual_review`: filtered for other reasons only.

Every row is in `*.all.tsv`, and also in the file named by its tier.

## Output columns

Rows are sorted by BAM header contig order, then position and strand; there
is one row per contig, position, and strand.

Shared artifact columns (in both modes): `support_mean_mapq local_mean_mapq
support_mean_nm local_mean_nm softclip_fraction secondary_endpoint_count
indel_position_0based indel_position_1based indel_type indel_length
indel_read_count indel_fraction known_indel_overlap`.

**Digenome** (one row per selected pair, keyed by the forward endpoint):
`sample analysis contig forward_position_0based forward_position_1based
reverse_position_0based reverse_position_1based tier filter_status
filter_reasons digenome_pair_score rgen_digenome_score forward_endpoint_count
forward_depth forward_fraction reverse_endpoint_count reverse_depth
reverse_fraction combined_endpoint_count combined_depth combined_fraction`,
the shared artifact columns, then `control_sample control_status
control_forward_endpoint_count control_forward_depth control_forward_fraction
control_reverse_endpoint_count control_reverse_depth control_reverse_fraction
control_combined_endpoint_count control_combined_depth
control_combined_fraction control_digenome_pair_score
control_rgen_digenome_score control_fold_enrichment control_fisher_p
control_fisher_q`.

**nDigenome** (one row per focal endpoint): `sample analysis contig
position_0based position_1based strand signal_classification tier
filter_status filter_reasons endpoint_count strand_depth endpoint_fraction
opposite_position_0based opposite_position_1based opposite_count
opposite_depth opposite_fraction`, the shared artifact columns, then
`control_sample control_status control_endpoint_count control_depth
control_fraction control_fold_enrichment control_fisher_p control_fisher_q`.

The BED name is `<sample>|<strand>|<class>` (`both` and `DSB` for Digenome),
the score is `min(1000, round(1000 × fraction))` using the combined or focal
fraction, and Digenome rows have strand `.`.

## Cutoff provenance

1. **Digenome reference settings.** The count, depth, and fraction defaults
   come from `digenome -G 0 -q 1 -f 5 -r 5 -d 10 -R 0.2 -s 2.5` in the
   [original Digenome distribution](http://www.rgenome.net/static/digenome-js/digenome).
   Its `-s 2.5` applies to the standalone tool's own score, which
   `rgen_digenome_score` reproduces (cross-checked against the v1.0
   executable). `digenome_pair_score` uses other math (strand fractions, no
   −1 terms, no five-position sum, divided by 4), so its default of 1.1 is
   2.5 moved to its scale: in a model of cut sites with 3–45 background reads
   per strand, balanced and 2:1 pairs, RGEN's 2.5 falls at a pair score of
   1.0–1.4, typically 1.1. The benchmark in `tests/benchmark/` reports how
   often the two cutoffs agree on real pairs.
2. **Published nDigenome focal thresholds.** At least 10 reads sharing a 5′
   endpoint and at least 20% local fraction come from Kim et al., *Unbiased
   investigation of specificities of prime editing systems in human cells*,
   Nucleic Acids Research (2020),
   [doi:10.1093/nar/gkaa764](https://doi.org/10.1093/nar/gkaa764).
3. **Pipeline safeguards.** The MAPQ floor, opposite-strand window, ambiguity
   thresholds, artifact window, clipping and indel limits, and control filters
   are implementation defaults, not taken from the nDigenome publication.

Calibrate production settings against known on-target sites, validated
off-target sites, matched negative controls, sequencing depth, library
preparation, and the nuclease used.
