# Digenome-seq and nDigenome-seq pipeline

Nextflow pipeline that calls nuclease cleavage sites from whole-genome
sequencing of digested genomic DNA: Illumina reads, as FASTQs or aligned BAMs,
or aligned ONT or PacBio long reads with `--long_reads`.

| Mode (`--analysis`) | Calls | Reads |
| --- | --- | --- |
| `digenome` (default) | forward and reverse endpoints paired into double-strand breaks | paired-end, single-end, or long reads |
| `ndigenome` | isolated strand endpoints: single-strand breaks or nicks | paired-end or long reads |

Steps: samplesheet check → bwa-mem2 index (shared cache) → fastp → bwa-mem2
alignment and duplicate marking → cleavage calling in parallel coordinate
chunks → per-sample finalize (pairing, sample-wide statistics, filters) →
MultiQC. Aligned BAMs skip trimming and alignment, are sorted if they need it,
and start at cleavage calling. How calls are made is in
[docs/cleavage_algorithm.md](docs/cleavage_algorithm.md); every parameter is
under [Parameters](#parameters).

## Requirements

- Nextflow 25.04 or later (Longleaf has 25.04.7) and Java 17 or later.
- On Longleaf: SLURM and Apptainer, plus the four images in this checkout's
  `containers/` directory (see [containers/README.md](containers/README.md)).
- Internet access on the launching node the first time, so Nextflow can
  download the pinned nf-schema plugin into `~/.nextflow/plugins`.

Check a Longleaf checkout:

```bash
command -v nextflow java apptainer
(cd containers && sha256sum -c checksums.sha256)
```

## Quick start

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --analysis digenome \
  --outdir results_digenome \
  -work-dir /work/groups/my_lab/digenome
```

`nextflow run /path/to/digenome-seq_pipeline --help` lists every parameter.
Unknown parameters and invalid values stop the run before anything is
submitted. Use separate output and work directories for separate runs; one
mode applies to every sample in a run.

## Samplesheet

```csv
sample,fastq_1,fastq_2,control,variant_vcf
Treated,/data/Treated_L001_R1.fastq.gz,/data/Treated_L001_R2.fastq.gz,Untreated,/data/donor.vcf.gz
Treated,/data/Treated_L002_R1.fastq.gz,/data/Treated_L002_R2.fastq.gz,Untreated,/data/donor.vcf.gz
Untreated,/data/Untreated_R1.fastq.gz,/data/Untreated_R2.fastq.gz,,
```

| Column | Required | Meaning |
| --- | --- | --- |
| `sample` | yes | Sample name: letters, digits, `.`, `_`, `-` |
| `fastq_1` | yes | R1 (or single-end) `.fastq.gz` / `.fq.gz` |
| `fastq_2` | yes | R2; blank only for single-end Digenome data |
| `control` | no | Name of another row's sample to use as the matched control |
| `variant_vcf` | no | bgzipped `.vcf.gz` with a `.tbi` or `.csi` index, for known-indel annotation |

- Repeat a sample name for more lanes; lanes are concatenated before trimming.
- A control is aligned but not called. It needs its own rows, with a blank
  `control`; one control can serve several samples.
- `control` and `variant_vcf` must be the same on every row of a sample.
- Each FASTQ file may appear only once (so R1 and R2 must differ).
- Unknown columns are rejected. Every problem in the sheet is reported at once.

A template is in `assets/samplesheet_template.csv`.

### Aligned BAMs

Reads you have already aligned can be given as BAMs instead: a `bam` column
replaces `fastq_1` and `fastq_2`. Trimming, alignment, and the reference cache
are skipped, and `--genome` isn't needed. One sheet lists either FASTQs or
BAMs, not both.

```csv
sample,bam,control,variant_vcf
Treated,/data/Treated.sorted.bam,Untreated,
Untreated,/data/Untreated.sorted.bam,,
```

- Each sample is one row and one aligned BAM; merge a sample's runs first.
- A BAM that is coordinate-sorted with `<bam>.bai` or `<bam>.csi` beside it
  is used as is. Any other is sorted and indexed first, into the work
  directory; your file is only read. Index your sorted BAMs to skip that step.
- Unaligned BAMs, such as raw basecaller output, are rejected.
- The BAMs are used as they are: reads flagged as duplicates are not counted,
  so mark duplicates first (for example with `samtools markdup`) if your
  aligner didn't. nDigenome still needs paired-end short reads.

## Long reads

`--long_reads` calls ONT or PacBio reads, given as [aligned BAMs](#aligned-bams)
(for example from minimap2).

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --long_reads \
  --input long_reads.csv \
  --analysis digenome \
  --outdir results_long_reads \
  -work-dir /work/groups/my_lab/long_reads
```

Both aligned ends of each primary read count, and the soft-clip and indel
limits become 1.0. [docs/cleavage_algorithm.md](docs/cleavage_algorithm.md#long-reads)
shows what cuts and nicks look like in long reads.

## Profiles and configuration

| Profile | Use |
| --- | --- |
| `longleaf` | SLURM, Apptainer, Longleaf binds and resource caps |
| `apptainer` | Local run with Apptainer |
| `slurm` | SLURM without containers |
| `test` | Tiny synthetic genome; list it last, e.g. `-profile apptainer,test` |

Resources and container images for each step are in `conf/base.config` and
apply in every profile. FASTP, ALIGN, and the two calling steps are retried
once with twice the memory and time if they are killed for running out of
either. Override resources, images, or genomes with your own file and
`-c custom.config`.

Configured genomes are `GRCh38`, `GRCh37`, and `GRCm39` (paths in
`nextflow.config`). `--genome` also takes their aliases, in any case: `hg38`
and `human_hg38`, `hg19` and `human_hg19`, `mm39` and `mouse_mm39`. All names
for a genome share one index. To add a genome:

```groovy
// custom.config
params.genomes = [MyGenome: [fasta: '/path/to/MyGenome.fa', aliases: ['mine']]]
```

then run with `-c custom.config --genome MyGenome` (or `--genome mine`).

## Parameters

Pipeline parameters take two hyphens (`--analysis ndigenome`); Nextflow's own
options take one (`-profile`). A boolean is turned on by its name alone
(`--long_reads`) or set with `true` or `false`. Unknown parameters and invalid
values stop the run at launch, and `--help` lists every parameter. The math
behind the cutoffs is in [docs/cleavage_algorithm.md](docs/cleavage_algorithm.md).

### Inputs and outputs

| Parameter | Default | Description |
| --- | --- | --- |
| `--input` | required | Samplesheet CSV of FASTQs or aligned BAMs; see [Samplesheet](#samplesheet). |
| `--genome` | none | Genome to align FASTQs to: a configured name or alias, in any case (`GRCh38` or `hg38`). Required for a FASTQ samplesheet; not used for BAMs. |
| `--analysis` | `digenome` | `digenome` pairs forward and reverse endpoints into double-strand breaks; `ndigenome` calls isolated strand endpoints (single-strand breaks, nicks). One mode applies to every sample. |
| `--outdir` | `results` | Where results are published. Use a separate one for each run. |
| `--genome_blacklist` | none | BED or BED.gz of regions to skip, in the BAMs' coordinates and contig names (see [Blacklist](#blacklist)). |
| `--ref_cache` | `reference_cache/` in the checkout | Shared bwa-mem2 index cache (see [Reference cache](#reference-cache)). Point it outside the checkout so a fresh clone reuses the index. |
| `--publish_trimmed_fastqs` | `false` | Also publish the fastp-trimmed FASTQs. |

### Run options

| Parameter | Default | Description |
| --- | --- | --- |
| `--long_reads` | `false` | Count both aligned ends of long ONT or PacBio reads, given as aligned BAMs (see [Long reads](#long-reads)). Also sets the soft-clip and indel limits to 1.0. |
| `--keep_multimappers` | `false` | Let MAPQ-0 primary alignments count, for repetitive regions (see [Multimappers](#multimappers)). |
| `--cleavage_chunks` | `8` | Coordinate chunks per sample, each called by a one-CPU task. Changes runtime, never results. |

### Trimming (FASTQ input only)

| Parameter | Default | Description |
| --- | --- | --- |
| `--fastp_qualified_quality_phred` | `20` | Minimum Phred quality for a base to count as qualified. |
| `--fastp_length_required` | `30` | Drop reads shorter than this after trimming. |
| `--fastp_low_complexity_filter` | `true` | Drop low-complexity reads in Digenome mode. Always off for nDigenome and with `--keep_multimappers`. |
| `--fastp_complexity_threshold` | `30` | Complexity threshold, in percent, for that filter. |
| `--fastp_extra_args` | empty | Extra arguments appended to the fastp command. |

Adapter detection is on for paired-end reads, and poly-G and poly-X trimming
for all reads.

### Digenome calling

A forward endpoint `f` and a reverse endpoint `r` can pair when
`|r - (f - overhang)| <= pair_window`. Each cutoff is strict: a value must be
greater than it, so a count cutoff of 5 needs at least 6 reads. Endpoints with
fewer reads than that are never reported.

| Parameter | Default | Description |
| --- | --- | --- |
| `--digenome_overhang` | `0` | Expected offset between the forward and reverse endpoints of one cut; the reverse endpoint is expected at `f - overhang`. |
| `--digenome_pair_window` | `2` | How many bases a reverse endpoint may sit from its expected position. |
| `--digenome_min_mapq` | `1` | Minimum MAPQ of counted reads (at least this). 0 with `--keep_multimappers`. |
| `--digenome_forward_cutoff` | `5` | Reads ending at the forward endpoint must be more than this. |
| `--digenome_reverse_cutoff` | `5` | Reads ending at the reverse endpoint must be more than this. |
| `--digenome_depth_cutoff` | `10` | Reads covering each endpoint must be more than this. |
| `--digenome_fraction_cutoff` | `0.2` | Each endpoint's share of the reads covering it must be more than this. |
| `--digenome_pair_score_cutoff` | `1.1` | `digenome_pair_score` (forward fraction × reverse fraction × the two counts summed ÷ 4) must be more than this. 1.1 is RGEN's cutoff of 2.5 on this score's scale; the RGEN-style score is reported for comparison and never filtered on. |

### nDigenome calling

Each threshold is inclusive: a value equal to it passes.

| Parameter | Default | Description |
| --- | --- | --- |
| `--ndigenome_min_count` | `10` | Reads ending at the endpoint must be at least this. |
| `--ndigenome_min_fraction` | `0.2` | Those reads' share of the same-strand reads covering the endpoint must be at least this. |
| `--ndigenome_min_mapq` | `1` | Minimum MAPQ of counted reads (at least this). 0 with `--keep_multimappers`. |
| `--ndigenome_opposite_window` | `5` | Bases on each side searched for an endpoint on the opposite strand. |
| `--ndigenome_ambiguous_min_count` | `3` | An opposite-strand endpoint with at least this many reads makes the call AMBIGUOUS. |
| `--ndigenome_ambiguous_min_fraction` | `0.05` | An opposite-strand endpoint with at least this share of its depth makes the call AMBIGUOUS. |

The strongest opposite-strand endpoint in the window sets the class:
POSSIBLE_DSB if it reaches `--ndigenome_min_count` and
`--ndigenome_min_fraction`, AMBIGUOUS if it reaches either ambiguous limit, and
SSB otherwise. Only unfiltered SSB rows are high confidence.

### Artifact and control filters

These apply in both modes. Each failed check adds its reason to the row's
`filter_reasons` and keeps the row out of the high-confidence calls.

| Parameter | Default | Description |
| --- | --- | --- |
| `--cleavage_artifact_window` | `10` | Bases on each side of an endpoint checked for MAPQ, mismatches, indels, clipping, and known indels. |
| `--cleavage_max_softclip_fraction` | `0.2` | HIGH_5P_SOFTCLIP when at least this share of the supporting reads are clipped at the end that forms the endpoint (the 5′ end for short reads). 1.0 with `--long_reads`. |
| `--cleavage_max_indel_fraction` | `0.2` | NEARBY_INDEL when at least this share of the reads in the window carry an insertion or deletion there. 1.0 with `--long_reads`. |
| `--cleavage_min_support_mean_mapq` | `10` | LOW_SUPPORT_MAPQ when the supporting reads' mean MAPQ is below this. 0 with `--keep_multimappers`. |
| `--cleavage_control_min_depth` | `1` | INSUFFICIENT_CONTROL_COVERAGE when fewer control reads than this cover the site. |
| `--cleavage_control_max_fraction` | `0.05` | HIGH_CONTROL_FRACTION when the control's share of reads ending at the site is above this. |
| `--cleavage_control_min_fold` | `5.0` | LOW_CONTROL_FOLD when the treated-over-control enrichment is below this. |
| `--cleavage_control_max_q` | `0.05` | CONTROL_Q_FAIL when the Fisher test's q-value, corrected across the whole sample, is above this. |

The four control filters apply only to samples that name a control.

### Multimappers

`--keep_multimappers` runs bwa-mem2 with `-a`, lets MAPQ-0 primary alignments
count (both minimum MAPQs and the support mean-MAPQ filter become 0), and turns
off fastp's low-complexity filter. Each read still counts once, at the
aligner's primary placement; secondary and supplementary alignments stay
diagnostic. Support can therefore be split across equivalent repeat copies.
With a samplesheet of BAMs only the MAPQ changes apply, since the reads come
aligned.

### Blacklist

`--genome_blacklist` takes 0-based half-open BED rows with the BAMs' contig
names; overlapping and adjacent rows are merged. Use a list made for the same
genome build as the BAMs: a list for another build either stops the run or
masks the wrong regions. Malformed rows, unknown contigs, out-of-range
coordinates, or an empty file stop the run. Endpoints inside the blacklist are
never called and never count as opposite-strand evidence. The QC JSON records
the file name, SHA-256, merged interval count, and excluded bases. The
pipeline never downloads or chooses a blacklist for you.

### Controls and variants

Controls are optional; rows without one are `UNCONTROLLED` and are not
filtered for that. Thresholds such as `--ndigenome_min_count` apply to each
treated sample on its own; a control is only compared at the same position.
With a control, each call gets fold enrichment, a Fisher exact p-value, and a
sample-wide Benjamini–Hochberg q-value; if control depth is below
`--cleavage_control_min_depth` the row is `INSUFFICIENT_CONTROL_COVERAGE` and
filtered. A VCF must declare every analyzed contig with the BAM's names and
lengths; a mismatch such as `chr1` versus `1` stops the run.

### Nextflow options

| Option | Description |
| --- | --- |
| `-profile` | Execution setup: `longleaf`, `apptainer`, `slurm`, or `test` (see [Profiles and configuration](#profiles-and-configuration)). |
| `-work-dir`, `-w` | Task work directory. Use a separate one for each run. |
| `-resume` | Reuse finished tasks from the same work directory. |
| `-c` | Add a configuration file, for custom genomes, resources, or images. |
| `-params-file` | Read parameters from a YAML or JSON file. |

The execution report, timeline, trace, and DAG are always written to
`<outdir>/pipeline_info/`.

### Removed in the 2026 rebuild

- `--max_memory`, `--max_cpus`, `--max_time`: set `process.resourceLimits` in a
  configuration file.
- `--containers`, `--container_bind_paths`: set images in `conf/base.config`
  and binds in a profile.
- `--index_lock_timeout_seconds`, `--index_stale_lock_seconds`: fixed at 48
  hours.
- `--publish_concat_fastqs`: lanes are concatenated inside FASTP.
- The samplesheet `lane` column: repeat the sample name on each lane's row.

## Outputs

```text
<outdir>/
├── <analysis>/
│   ├── <sample>.<analysis>.all.tsv              every evaluated row
│   ├── <sample>.<analysis>.high_confidence.tsv  rows that pass every filter
│   ├── <sample>.<analysis>.bed                  the same calls as BED6
│   ├── <sample>.<analysis>.artifact.tsv         filtered, with artifact evidence
│   ├── <sample>.<analysis>.manual_review.tsv    filtered for other reasons
│   ├── <sample>.<analysis>.qc.json
│   └── <sample>.<analysis>_mqc.tsv
├── bam/  fastp/  qc/                            short reads only
├── multiqc/
├── trimmed_fastqs/                              only with --publish_trimmed_fastqs
└── pipeline_info/
    ├── analysis_parameters.json                 resolved settings, written at launch
    ├── container_checksums.sha256, container_sources.tsv
    ├── cleavage_chunks/<sample>.cleavage_chunks.tsv
    └── Nextflow report, timeline, trace, and DAG
```

The `tier` column names the file each row is in. Digenome and nDigenome tables
have different columns, and rows are sorted by the BAM's contig order; both
are listed in [docs/cleavage_algorithm.md](docs/cleavage_algorithm.md#output-columns).
`digenome_pair_score` drives Digenome pairing and filtering;
`rgen_digenome_score` reproduces the standalone CRISPR RGEN Tools v1.0 score
for comparison with historical results only.

## Reference cache

bwa-mem2 indexes are built once and shared:

```text
<ref_cache>/<genome>/<fasta_sha256>/<bwa_mem2_version>/bwamem2/
```

`--ref_cache` defaults to `reference_cache/` in the checkout being run. An
index is reused only when the genome name, the FASTA's SHA-256, the bwa-mem2
version, and every index file match. Concurrent runs share a lock: one builds
while the others wait, then reuse it. Permissions are set explicitly so
several users can share the cache whatever their umask:

```text
<ref_cache>/                         1777
<ref_cache>/<genome>/                1777
<ref_cache>/<genome>/<fasta_sha256>/ 1777
.../<bwa_mem2_version>/              1777
.../<bwa_mem2_version>/bwamem2/      755
.../bwamem2/index files              644
```

The sticky bit lets users add cache entries without deleting or renaming each
other's. Lock directories and their owner records are world-readable so a
waiting user can see who is building. An extra setgid bit (modes such as 3777
or 2755) is accepted.

## Testing and validation

`./tests/run_tests.sh` runs the local suite; [docs/testing.md](docs/testing.md)
explains it, the smoke test, and what only a Longleaf run can show. Before
production use, run both modes through the Longleaf smoke test and validate
with known-positive material, untreated controls, and representative library
preparations, plus a repetitive-region truth set for `--keep_multimappers`.
