# Digenome-seq and nDigenome-seq pipeline

Nextflow pipeline that calls nuclease cleavage sites from whole-genome
sequencing of digested genomic DNA: Illumina reads, or aligned ONT or PacBio
long reads with `--long_reads`.

| `--analysis` | Calls | Reads |
| --- | --- | --- |
| `digenome` (default) | forward and reverse endpoints paired into double-strand breaks | paired-end, single-end, or long reads |
| `ndigenome` | isolated strand endpoints: single-strand breaks or nicks | paired-end or long reads |

Steps: samplesheet check → bwa-mem2 index (shared cache) → fastp → bwa-mem2
alignment and duplicate marking → cleavage calling in parallel coordinate
chunks → per-sample finalize (pairing, sample-wide statistics, filters) →
MultiQC. Long reads arrive aligned, are sorted if they need it, and start at
cleavage calling. How calls are made is in
[docs/cleavage_algorithm.md](docs/cleavage_algorithm.md); every parameter is in
[docs/parameters.md](docs/parameters.md).

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

## Long reads

`--long_reads` calls ONT or PacBio reads from BAMs you have already aligned,
for example with minimap2. Trimming, alignment, and the reference cache are
skipped, and `--genome` isn't needed.

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --long_reads \
  --input long_reads.csv \
  --analysis digenome \
  --outdir results_long_reads \
  -work-dir /work/groups/my_lab/long_reads
```

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
- `control` and `variant_vcf` work as they do for short reads.
- Reads flagged as duplicates are not counted, as with short reads; marking
  them is up to you.
- Both aligned ends of each primary read count, and the soft-clip and indel
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

## Key options

| Option | Default | Effect |
| --- | --- | --- |
| `--long_reads` | off | Call aligned long-read BAMs (see [Long reads](#long-reads)) |
| `--keep_multimappers` | off | Count MAPQ-0 primary alignments for repetitive regions (see below) |
| `--cleavage_chunks` | 8 | Chunks per sample, each a one-CPU task; changes runtime, never results |
| `--genome_blacklist` | none | BED/BED.gz of regions to skip (see below) |
| `--publish_trimmed_fastqs` | off | Also publish the trimmed FASTQs |

**Multimappers.** `--keep_multimappers` runs bwa-mem2 with `-a`, lets MAPQ-0
primary alignments count (both minimum MAPQs and the support mean-MAPQ filter
become 0), and turns off fastp's low-complexity filter. Each read still counts
once, at BWA's primary placement; secondary and supplementary alignments stay
diagnostic. Support can therefore be split across equivalent repeat copies.
With `--long_reads` only the MAPQ changes apply, since the BAMs come
aligned.

**Blacklist.** `--genome_blacklist` takes 0-based half-open BED rows with the
BAM's contig names; overlapping and adjacent rows are merged. Malformed rows,
unknown contigs, out-of-range coordinates, or an empty file stop the run.
Endpoints inside the blacklist are never called and never count as
opposite-strand evidence. The QC JSON records the file name, SHA-256, merged
interval count, and excluded bases. The pipeline never downloads or chooses a
blacklist for you.

**Controls and variants.** Controls are optional; rows without one are
`UNCONTROLLED` and are not filtered for that. With a control, each call gets
fold enrichment, a Fisher exact p-value, and a sample-wide
Benjamini–Hochberg q-value; if control depth is below
`--cleavage_control_min_depth` the row is `INSUFFICIENT_CONTROL_COVERAGE` and
filtered. A VCF must declare every analyzed contig with the BAM's names and
lengths; a mismatch such as `chr1` versus `1` stops the run.

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
