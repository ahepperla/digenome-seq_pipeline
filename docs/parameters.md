# Parameter reference

This page explains every pipeline parameter. Defaults live in
`nextflow.config`; types, ranges, and short descriptions live in
`nextflow_schema.json` (`--help` prints them). Pipeline parameters use two leading hyphens and underscores:

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --analysis ndigenome \
  --cleavage_chunks 64
```

Boolean parameters can be enabled by name, such as `--keep_multimappers`, or
set explicitly with `true` or `false`. Quote paths and free-form strings when
they contain spaces.

`nextflow_schema.json` defines every accepted pipeline parameter, type,
default, required value, choice, and numeric range. The nf-schema plugin
(pinned 2.5.1) validates parameters at launch: it rejects unknown parameters
and invalid values, and checks that required parameters are supplied.

## Inputs and outputs

| Parameter | Type | Default | Required | Description |
| --- | --- | --- | --- | --- |
| `--input` | Path | `null` | Yes | Samplesheet CSV: `sample`, `fastq_1`, `fastq_2`, and optional `control` and `variant_vcf`; or, for aligned BAMs (required with `--long_reads`), `sample`, `bam`, and optional `control` and `variant_vcf`. See README.md for the rules. |
| `--genome` | String | `null` | For a samplesheet of FASTQs | Name or alias of a genome in params.genomes to align against, in any case (for example GRCh38 or hg38). A samplesheet of BAMs is already aligned, so it isn't needed then. |
| `--analysis` | Choice | `digenome` | No | Calling mode: `digenome` pairs forward and reverse endpoints into DSBs; `ndigenome` calls isolated strand endpoints (SSBs, nicks). One mode applies to the complete run. |
| `--outdir` | Path | `results` | No | Published results directory. Use a separate one for each run. |
| `--genome_blacklist` | Path | `null` | No | Optional BED or BED.gz of regions to skip, with BAM contig names and 0-based half-open coordinates. |
| `--ref_cache` | Path | `${projectDir}/reference_cache` | No | Shared bwa-mem2 index cache, keyed by genome, FASTA SHA-256, and bwa-mem2 version. |
| `--publish_trimmed_fastqs` | Boolean | `false` | No | Also publish the fastp-trimmed FASTQs. |

`params.genomes` maps genome names to a FASTA path and optional aliases. Add or
override entries in a Nextflow configuration file rather than on the command
line. Repository defaults are `GRCh38` (aliases `hg38`, `human_hg38`),
`GRCh37` (`hg19`, `human_hg19`), and `GRCm39` (`mm39`, `mouse_mm39`). The index
cache and the run record use the name, so a genome's aliases share one index.

## Execution

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `--keep_multimappers` | Boolean | `false` | Run bwa-mem2 with `-a`, count MAPQ-0 primary alignments (both minimum MAPQs and the support mean-MAPQ filter become 0), and turn off fastp low-complexity filtering. Each read still counts once, at its primary placement. With a samplesheet of BAMs only the MAPQ changes apply. |
| `--cleavage_chunks` | Integer | `8` | Coordinate chunks per sample, each called by a one-CPU task. Changes runtime, never results. |
| `--long_reads` | Boolean | `false` | Call long-read (ONT/PacBio) samples from aligned BAMs given in the samplesheet's `bam` column; trimming and alignment are skipped, and BAMs that aren't coordinate-sorted and indexed are sorted first. Both aligned ends of each read count (see [cleavage_algorithm.md](cleavage_algorithm.md#long-reads)), and `--cleavage_max_softclip_fraction` and `--cleavage_max_indel_fraction` become 1.0 because nanopore reads routinely carry small indels and clipped ends. |

## Read trimming (fastp)

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `--fastp_qualified_quality_phred` | Integer | `20` | Minimum Phred score for a qualified base. |
| `--fastp_length_required` | Integer | `30` | Discard reads shorter than this after trimming. |
| `--fastp_low_complexity_filter` | Boolean | `true` | Low-complexity filtering in Digenome mode. Always off in nDigenome mode and with --keep_multimappers. |
| `--fastp_complexity_threshold` | Integer | `30` | fastp complexity threshold (percent) when low-complexity filtering is on. |
| `--fastp_extra_args` | String | Empty | Extra arguments appended to the fastp command. |

Adapter detection is enabled for paired-end reads. Poly-G and poly-X trimming
are enabled for both paired-end and single-end reads.

## Digenome calling

These parameters affect `--analysis digenome`. Forward and reverse endpoints
are eligible to pair when:

```text
abs(reverse_position - (forward_position - overhang)) <= pair_window
```

The filtering score is:

```text
digenome_pair_score = forward_fraction * reverse_fraction
                    * (forward_endpoint_count + reverse_endpoint_count) / 4
```

| Parameter | Type | Default | Comparison | Description |
| --- | --- | --- | --- | --- |
| `--digenome_overhang` | Integer | `0` | Coordinate adjustment | Expected forward-minus-reverse endpoint offset; the expected reverse endpoint is forward - overhang. |
| `--digenome_pair_window` | Integer | `2` | Within `+/-` window | A forward and reverse endpoint can pair when `\|reverse - (forward - overhang)\| <= this`. |
| `--digenome_min_mapq` | Integer | `1` | Alignment MAPQ `>=` value | Minimum MAPQ (>=) of counted alignments. 0 with --keep_multimappers. |
| `--digenome_forward_cutoff` | Integer | `5` | Forward count `>` value | Forward endpoint count must be > this. |
| `--digenome_reverse_cutoff` | Integer | `5` | Reverse count `>` value | Reverse endpoint count must be > this. |
| `--digenome_depth_cutoff` | Integer | `10` | Each strand depth `>` value | Each strand's depth must be > this. |
| `--digenome_fraction_cutoff` | Number | `0.2` | Each endpoint fraction `>` value | Each strand's endpoint fraction must be > this. |
| `--digenome_pair_score_cutoff` | Number | `1.1` | Pair score `>` value | digenome_pair_score must be > this. 1.1 is the standalone RGEN cutoff of 2.5 moved to this score's scale (see below); the two scores' cutoffs are not interchangeable. |

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

## nDigenome calling

These parameters affect `--analysis ndigenome`.

| Parameter | Type | Default | Comparison | Description |
| --- | --- | --- | --- | --- |
| `--ndigenome_min_count` | Integer | `10` | Focal count `>=` value | Focal endpoint count must be >= this. An opposite endpoint reaching this and the minimum fraction makes the call POSSIBLE_DSB. |
| `--ndigenome_min_fraction` | Number | `0.2` | Focal fraction `>=` value | Focal endpoint fraction of same-strand depth must be >= this. |
| `--ndigenome_min_mapq` | Integer | `1` | Alignment MAPQ `>=` value | Minimum MAPQ (>=) of counted alignments. 0 with --keep_multimappers. |
| `--ndigenome_opposite_window` | Integer | `5` | `+/-` window | Search this many bases on each side for opposite-strand endpoints. |
| `--ndigenome_ambiguous_min_count` | Integer | `3` | Opposite count `>=` value | An opposite endpoint with count >= this makes the call AMBIGUOUS. |
| `--ndigenome_ambiguous_min_fraction` | Number | `0.05` | Opposite fraction `>=` value | An opposite endpoint with fraction >= this makes the call AMBIGUOUS. |

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
| `--cleavage_artifact_window` | Integer | `10` | Measurement window | Bases on each side of an endpoint used for MAPQ, mismatch, indel, clipping, and known-indel checks. |
| `--cleavage_max_softclip_fraction` | Number | `0.2` | Soft-clipped fraction `>=` value | HIGH_5P_SOFTCLIP when the fraction of supporting reads clipped at the end that forms the endpoint (the 5' end for short reads) is >= this. |
| `--cleavage_max_indel_fraction` | Number | `0.2` | Local indel fraction `>=` value | NEARBY_INDEL when the fraction of local alignments with a nearby indel is >= this. |
| `--cleavage_min_support_mean_mapq` | Number | `10` | Mean MAPQ `<` value | LOW_SUPPORT_MAPQ when supporting reads' mean MAPQ is < this. 0 with --keep_multimappers. |
| `--cleavage_control_min_depth` | Integer | `1` | Control depth `<` value | INSUFFICIENT_CONTROL_COVERAGE when the control depth is < this. |
| `--cleavage_control_max_fraction` | Number | `0.05` | Control fraction `>` value | HIGH_CONTROL_FRACTION when the control endpoint fraction is > this. Equality passes. |
| `--cleavage_control_min_fold` | Number | `5.0` | Fold enrichment `<` value | LOW_CONTROL_FOLD when the treated/control fold enrichment is < this. Equality passes. |
| `--cleavage_control_max_q` | Number | `0.05` | Fisher q-value `>` value | CONTROL_Q_FAIL when the sample-wide Benjamini-Hochberg q-value is > this. Equality passes. |

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

## Removed in the 2026 rebuild

These parameters were used in the previous pipeline and are no longer available.
Resource limits, container image paths, and index cache timeouts are now
configured through configuration files instead of pipeline parameters:

- `max_memory`: use `resourceLimits` in a configuration file
- `max_cpus`: use `resourceLimits` in a configuration file
- `max_time`: use `resourceLimits` in a configuration file
- `containers`: set process container paths in `conf/base.config` or a custom config
- `container_bind_paths`: the `longleaf` profile binds required paths; customize in a config
- `index_lock_timeout_seconds`: fixed at 48 hours in `bin/prepare_bwamem2_index.sh`
- `index_stale_lock_seconds`: fixed at 48 hours in `bin/prepare_bwamem2_index.sh`
- `publish_concat_fastqs`: lane concatenation happens inside the FASTP process
- the samplesheet `lane` column: repeat the sample name on each lane's row instead

## Nextflow runtime options

These are Nextflow options, not pipeline parameters. They use one leading
hyphen and are interpreted by Nextflow itself.

| Option | Default | Description |
| --- | --- | --- |
| `-profile` | None | Selects execution configuration: `longleaf`, `apptainer`, `slurm`, or `test` (list `test` last, e.g. `apptainer,test`). Production runs should use `longleaf`. |
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
  --genome GRCh38 \
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
  --genome GRCh38 \
  --analysis digenome \
  --digenome_overhang 4 \
  --digenome_pair_window 2 \
  --outdir results_digenome \
  -work-dir /work/groups/example/digenome
```
