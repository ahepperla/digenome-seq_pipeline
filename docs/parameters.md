# Parameter reference

This page is the authoritative lookup for pipeline parameters in
`nextflow.config`. Pipeline parameters use two leading hyphens and underscores:

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --input samplesheet.csv \
  --genome hg38 \
  --analysis ndigenome \
  --cleavage_chunks 64
```

Boolean parameters can be enabled by name, such as `--keep_multimappers`, or
set explicitly with `true` or `false`. Quote paths and free-form strings when
they contain spaces.

`nextflow_schema.json` defines every accepted pipeline parameter, type,
default, required value, choice, and numeric range. A preflight task rejects
unknown parameters and invalid values, checks configured input/reference/
blacklist/container/bind paths, and checks output/cache accessibility before
samplesheet validation or reference-index preparation can begin. Its report
is published as `<outdir>/pipeline_info/preflight.ready.json`.

## Required inputs and references

| Parameter | Type | Default | Required | Description |
| --- | --- | --- | --- | --- |
| `--input` | Path | `null` | Yes | Input samplesheet CSV. See the samplesheet section of the main README for columns and validation rules. |
| `--genome` | String | `null` | Yes | Configured genome name or alias. Defaults supplied by the repository are `GRCh38`/`hg38`, `GRCh37`/`hg19`, and `GRCm39`/`mm39`. Matching is case-insensitive. |
| `--analysis` | Choice | `digenome` | No | Calling mode: `digenome` for paired DSB endpoints or `ndigenome` for strand-specific SSB endpoints. One mode applies to the complete run. |
| `--outdir` | Path | `results` | No | Published output directory. Use a separate output directory for each independent run. |
| `--genome_blacklist` | Path | `null` | No | Optional BED or BED.gz mask using BAM-matching contig names and 0-based half-open coordinates. Masked regions are skipped during candidate scanning. Masked focal endpoints, Digenome pair endpoints, and nDigenome opposite-strand evidence are excluded. |
| `--ref_cache` | Path | `${projectDir}/reference_cache` | No | Shared, content-addressed bwa-mem2 index cache stored under `<ref_cache>/<genome>/<fasta_sha256>/<bwa_mem2_version>/bwamem2/`. The helper applies explicit shared permissions independent of user umask. The Longleaf profile also uses this project-relative default. |
| `params.genomes` | Map | Repository genome map | Configuration | Maps genome names to a FASTA path and optional aliases. Add or override entries in a Nextflow configuration file rather than on the command line. |
| `params.containers` | Map | `${projectDir}/containers/*.sif` | Configuration | Maps the `python`, `fastp`, `align`, `cleavage`, and `multiqc` process groups to container images. Every entry is required. |
| `params.container_bind_paths` | List | `[]` | Configuration | Host paths bound into Apptainer containers. The Longleaf profile sets `/proj`, `/work`, `/users`, `/overflow`, and `/nas`. |
| `--index_lock_timeout_seconds` | Integer seconds | `172800` | No | Maximum time to wait for another process that is building the same bwa-mem2 index. The default is 48 hours. |
| `--index_stale_lock_seconds` | Integer seconds | `172800` | No | Age at which an abandoned index-build lock becomes eligible for recovery. The default is 48 hours. |

## Execution and publishing

| Parameter | Type | Default | Applies to | Description |
| --- | --- | --- | --- | --- |
| `--max_memory` | Memory | `3041.GB` | SLURM profiles | Executor-wide resource cap. Individual process requests remain defined in `conf/base.config`. |
| `--max_cpus` | Integer | `256` | SLURM profiles | Executor-wide CPU cap, not the CPU request for an individual task. |
| `--max_time` | Duration | `240.h` | SLURM profiles | Executor-wide wall-time cap. Individual process limits remain defined in `conf/base.config`. |
| `--keep_multimappers` | Boolean | `false` | Both modes | Runs bwa-mem2 with `-a`, permits MAPQ 0 primary alignments, sets both caller minimum MAPQ values to 0, disables the support mean-MAPQ filter, and disables fastp low-complexity filtering. Secondary and supplementary alignments remain diagnostic only and are not counted as independent cleavage support. |
| `--cleavage_chunks` | Integer | `8` | Both modes | Requested coordinate chunks per sample and maximum concurrent cleavage-caller tasks. Each caller uses one CPU. Work estimates exclude supplied blacklist bases, so masked spans do not consume chunk capacity. The planner can emit fewer nonempty chunks when the callable data cannot support the requested count. |
| `--publish_concat_fastqs` | Boolean | `false` | Both modes | Copies lane-concatenated FASTQs to `<outdir>/concat_fastqs`. Concatenated files are always created in work storage for downstream processing. |
| `--publish_trimmed_fastqs` | Boolean | `false` | Both modes | Copies fastp-trimmed FASTQs to `<outdir>/trimmed_fastqs`. Trimmed files are always created in work storage for alignment. |

Chunk ownership is complete, gap-free, and nonoverlapping. Each chunk scans a
padded region so nearby strand evidence is preserved, but each call has one
owner. Final filtering and Benjamini-Hochberg q-values are calculated across
the complete sample after all chunks are merged. Blacklisted spans remain in
the ownership map for continuity, but contribute zero estimated callable work.

## fastp preprocessing

| Parameter | Type | Default | Applies to | Description |
| --- | --- | --- | --- | --- |
| `--fastp_qualified_quality_phred` | Integer | `20` | Both modes | Minimum Phred score used by fastp to classify a base as qualified. |
| `--fastp_length_required` | Integer bases | `30` | Both modes | Discards reads shorter than this length after trimming. |
| `--fastp_low_complexity_filter` | Boolean | `true` | Digenome only | Enables fastp low-complexity filtering in Digenome mode. It is automatically disabled in nDigenome mode and whenever `--keep_multimappers` is enabled. |
| `--fastp_complexity_threshold` | Integer percent | `30` | Digenome only | fastp complexity threshold used when low-complexity filtering is active. |
| `--fastp_extra_args` | String | Empty | Both modes | Additional arguments appended to the fastp command. This is an advanced escape hatch; avoid duplicating options already controlled by dedicated parameters. |

Adapter detection is enabled for paired-end reads. Poly-G and poly-X trimming
are enabled for both paired-end and single-end reads.

## Digenome parameters

These parameters affect `--analysis digenome`. Forward and reverse endpoints
are eligible to pair when:

```text
abs(reverse_position - (forward_position - overhang)) <= pair_window
```

The filtering score is:

```text
forward_fraction * reverse_fraction
* (forward_endpoint_count + reverse_endpoint_count) / 4
```

| Parameter | Type | Default | Pass condition | Description |
| --- | --- | --- | --- | --- |
| `--digenome_overhang` | Signed integer bases | `0` | Coordinate adjustment | Expected forward-minus-reverse endpoint offset. A positive value places the expected reverse endpoint at `forward - overhang`. |
| `--digenome_pair_window` | Nonnegative integer bases | `2` | Within `+/-` window | Positional tolerance around the expected reverse endpoint. For example, overhang `4` and window `2` accept forward-minus-reverse offsets from 2 through 6 bases. |
| `--digenome_min_mapq` | Integer | `1` | Alignment MAPQ `>=` value | Minimum MAPQ for primary alignments contributing to Digenome candidates and metrics. Automatically becomes `0` with `--keep_multimappers`. |
| `--digenome_forward_cutoff` | Integer count | `5` | Forward count `>` value | A cutoff of 5 requires at least 6 forward endpoint reads. |
| `--digenome_reverse_cutoff` | Integer count | `5` | Reverse count `>` value | A cutoff of 5 requires at least 6 reverse endpoint reads. |
| `--digenome_depth_cutoff` | Integer count | `10` | Each strand depth `>` value | Both forward and reverse local strand depths must exceed this cutoff. |
| `--digenome_fraction_cutoff` | Fraction | `0.20` | Each endpoint fraction `>` value | Both strand-specific endpoint fractions must exceed this cutoff. |
| `--digenome_pair_score_cutoff` | Number | `2.5` | Pair score `>` value | Minimum `digenome_pair_score`. |

Candidate pairs are selected with deterministic one-to-one matching. Pairs
that fail caller thresholds remain in the complete audit TSV with their
specific filter reasons.

The call TSV reports two Digenome scores:

- `digenome_pair_score` is the strand-specific score above and is the only
  score used for pairing priority and filtering.
- `rgen_digenome_score` reproduces the standalone CRISPR RGEN Tools v1.0
  five-position, total-depth calculation for comparison with historical
  results. It is not used for filtering.

Matched-control rows report the independently measured
`control_digenome_pair_score` and `control_rgen_digenome_score`.

## nDigenome parameters

These parameters affect `--analysis ndigenome`.

| Parameter | Type | Default | Comparison | Description |
| --- | --- | --- | --- | --- |
| `--ndigenome_min_count` | Integer count | `10` | Focal count `>=` value | Minimum number of same-strand reads ending at the focal coordinate. |
| `--ndigenome_min_fraction` | Fraction | `0.20` | Focal fraction `>=` value | Minimum fraction of local same-strand depth ending at the focal coordinate. |
| `--ndigenome_min_mapq` | Integer | `1` | Alignment MAPQ `>=` value | Minimum MAPQ for primary alignments contributing to nDigenome candidates and metrics. Automatically becomes `0` with `--keep_multimappers`. |
| `--ndigenome_opposite_window` | Nonnegative integer bases | `5` | `+/-` window | Searches for opposite-strand endpoints within this distance of the focal endpoint. |
| `--ndigenome_ambiguous_min_count` | Nonnegative integer count | `3` | Opposite count `>=` value | Weak opposite-strand count threshold. Meeting this threshold alone is enough to classify the row as `AMBIGUOUS` when it is not already `POSSIBLE_DSB`. |
| `--ndigenome_ambiguous_min_fraction` | Fraction | `0.05` | Opposite fraction `>=` value | Weak opposite-strand fraction threshold. Meeting this threshold alone is enough to classify the row as `AMBIGUOUS` when it is not already `POSSIBLE_DSB`. |

Opposite-strand classification uses:

```text
POSSIBLE_DSB:
    opposite_count >= ndigenome_min_count
    AND opposite_fraction >= ndigenome_min_fraction

AMBIGUOUS:
    opposite_count >= ndigenome_ambiguous_min_count
    OR opposite_fraction >= ndigenome_ambiguous_min_fraction

SSB:
    neither opposite-strand condition is met
```

Only unfiltered `SSB` rows enter the high-confidence nDigenome output.
`POSSIBLE_DSB` and `AMBIGUOUS` rows remain in the complete audit TSV.

## Shared artifact and control filters

| Parameter | Type | Default | Filter condition | Description |
| --- | --- | --- | --- | --- |
| `--cleavage_artifact_window` | Nonnegative integer bases | `10` | Measurement window | Radius around each endpoint used for local MAPQ, mismatch, indel, and soft-clipping measurements. It also contributes to automatic chunk padding. |
| `--cleavage_max_softclip_fraction` | Fraction | `0.20` | Soft-clipped fraction `>=` value | Adds `HIGH_5P_SOFTCLIP` when this fraction of supporting endpoint reads has a 5-prime soft clip. |
| `--cleavage_max_indel_fraction` | Fraction | `0.20` | Local indel fraction `>=` value | Adds `NEARBY_INDEL` when the fraction of local primary alignments with a nearby CIGAR indel reaches this value. |
| `--cleavage_min_support_mean_mapq` | Number | `10` | Mean MAPQ `<` value | Adds `LOW_SUPPORT_MAPQ` when the mean MAPQ of supporting endpoint reads is below this value. Automatically becomes `0` with `--keep_multimappers`. |
| `--cleavage_control_min_depth` | Integer count | `1` | Control depth `<` value | Adds `INSUFFICIENT_CONTROL_COVERAGE` when matched-control strand depth is below this value. It has no effect when no matched control is supplied. |
| `--cleavage_control_max_fraction` | Fraction | `0.05` | Control fraction `>` value | Adds `HIGH_CONTROL_FRACTION`. Equality passes. |
| `--cleavage_control_min_fold` | Number | `5.0` | Fold enrichment `<` value | Adds `LOW_CONTROL_FOLD` when the pseudocount-adjusted treated endpoint rate is not sufficiently enriched over control. Equality passes. |
| `--cleavage_control_max_q` | Fraction | `0.05` | Fisher q-value `>` value | Adds `CONTROL_Q_FAIL` when the sample-wide Benjamini-Hochberg-adjusted matched-control Fisher value exceeds this threshold. Equality passes. |

For matched controls, each failed condition is reported independently:

```text
HIGH_CONTROL_FRACTION:
    control_fraction > cleavage_control_max_fraction

LOW_CONTROL_FOLD:
    treated/control fold < cleavage_control_min_fold

CONTROL_Q_FAIL:
    control_fisher_q > cleavage_control_max_q
```

The fold enrichment uses:

```text
treated_rate = (treated_endpoint_count + 0.5) / (treated_depth + 1)
control_rate = (control_endpoint_count + 0.5) / (control_depth + 1)
fold = treated_rate / control_rate
```

The Fisher exact test itself uses the raw endpoint and non-endpoint counts.
Q-values are calculated after merging every chunk so multiple-testing
correction remains sample-wide.

## Nextflow runtime options

These are Nextflow options, not pipeline parameters. They use one leading
hyphen and are interpreted by Nextflow itself.

| Option | Default | Description |
| --- | --- | --- |
| `-profile` | None | Selects execution configuration. Available profiles are `standard`, `slurm`, `apptainer`, and `longleaf`; production runs should select the appropriate profile explicitly. |
| `-resume` | Off | Reuses compatible cached tasks from the selected work directory. Changed scripts, parameters, or inputs cause affected tasks to run again. |
| `-work-dir` | `work` | Nextflow task work directory. Use stable, separate locations for independent production runs. |
| `-c` | None | Adds a Nextflow configuration file. Useful for custom genomes, containers, resources, or site settings. |
| `-with-report` | Profile configured | Enables or overrides the execution report. The bundled configuration writes one under `<outdir>/pipeline_info`. |
| `-with-trace` | Profile configured | Enables or overrides the task trace. The bundled configuration writes one under `<outdir>/pipeline_info`. |
| `-with-timeline` | Profile configured | Enables or overrides the timeline. The bundled configuration writes one under `<outdir>/pipeline_info`. |
| `-with-dag` | Profile configured | Enables or overrides the workflow DAG. The bundled configuration writes one under `<outdir>/pipeline_info`. |

## Examples

High-parallelism nDigenome run with repetitive-region support and a supplied
blacklist:

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  -resume \
  --input samplesheet.csv \
  --genome hg38 \
  --analysis ndigenome \
  --keep_multimappers \
  --genome_blacklist /path/to/genome_blacklist.bed.gz \
  --cleavage_chunks 64 \
  --outdir results_ndigenome \
  -work-dir /work/groups/example/ndigenome
```

Digenome run expecting a four-base endpoint offset with two-base tolerance:

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --input samplesheet.csv \
  --genome hg38 \
  --analysis digenome \
  --digenome_overhang 4 \
  --digenome_pair_window 2 \
  --outdir results_digenome \
  -work-dir /work/groups/example/digenome
```
