# Digenome-seq and nDigenome-seq pipeline

Nextflow DSL2 pipeline for detecting nuclease cleavage sites from Illumina
whole-genome sequencing.

| Mode | Detects | Read requirement |
| --- | --- | --- |
| `digenome` | Paired strand endpoints representing double-strand breaks | Paired-end or single-end |
| `ndigenome` | Isolated strand endpoints representing single-strand breaks or nicks | Paired-end |

Both modes use the same alignment and cleavage-calling engine. They share
matched controls, known-indel annotation, artifact filters, parallel cleavage
calling, QC, and output formats.

## Requirements

- Nextflow >= 25.04 (Longleaf has 25.04.7)
- Java 17+
- SLURM and Apptainer on Longleaf only

The pipeline downloads the nf-schema plugin (pinned 2.5.1) on first run into
`~/.nextflow/plugins`, so the launching node needs internet once.

Container images provide Python 3.12, pysam 0.23.3, bwa-mem2, samtools, fastp,
and MultiQC. No Python module or virtual environment is required on the host.

## Installation and verification

SIF files are excluded from Git because of their size. Each complete installation
stores the validated images in:

```text
<pipeline-directory>/containers
```

Nextflow expects the files in the `containers/` directory of whichever pipeline
checkout is being run. A separate checkout must have these images copied or
provisioned there before running the pipeline.

Verify the required commands and the default GRCh38 reference:

```bash
command -v nextflow java apptainer
test -r /proj/seq/data/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa
```

Verify the shared Longleaf containers:

```bash
cd /path/to/digenome-seq_pipeline/containers
sha256sum -c checksums.sha256
```

Build the cleavage image only when it is missing or intentionally being replaced:

```bash
cd /path/to/digenome-seq_pipeline
./containers/build_cleavage.sh
```

## Quick start

Run Digenome-seq:

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --analysis digenome \
  --outdir results_digenome \
  -work-dir /work/groups/barc_scr/my_project/digenome
```

Run nDigenome-seq with repetitive-region support and 16 cleavage callers:

```bash
nextflow run /path/to/digenome-seq_pipeline \
  -profile longleaf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --analysis ndigenome \
  --keep_multimappers \
  --cleavage_chunks 16 \
  --outdir results_ndigenome \
  -work-dir /work/groups/barc_scr/my_project/ndigenome
```

Always use separate output and work directories for separate runs. One
analysis mode applies to every analyzed sample in that run. If a project
contains both DSB and SSB samples, use separate mode-specific samplesheets or
otherwise ensure each run contains only samples appropriate for that mode.
A shared control can appear in both samplesheets.

## Samplesheet

Minimal paired-end samplesheet:

```csv
sample,fastq_1,fastq_2
SampleA,/data/SampleA_R1.fastq.gz,/data/SampleA_R2.fastq.gz
```

With a matched control and optional known variants:

```csv
sample,fastq_1,fastq_2,control,variant_vcf
Treated,/data/Treated_R1.fastq.gz,/data/Treated_R2.fastq.gz,Untreated,/data/donor.vcf.gz
Untreated,/data/Untreated_R1.fastq.gz,/data/Untreated_R2.fastq.gz,,
```

Each row is one FASTQ pair. Repeat the same `sample` name to combine multiple
lanes or files before trimming and alignment:

```csv
sample,fastq_1,fastq_2,control
Treated,/data/Treated_L001_R1.fastq.gz,/data/Treated_L001_R2.fastq.gz,Untreated
Treated,/data/Treated_L002_R1.fastq.gz,/data/Treated_L002_R2.fastq.gz,Untreated
Untreated,/data/Untreated_R1.fastq.gz,/data/Untreated_R2.fastq.gz,
```

Columns:

| Column | Required | Meaning |
| --- | --- | --- |
| `sample` | Yes | Sample identifier used in filenames and reports |
| `fastq_1` | Yes | R1 or single-end FASTQ |
| `fastq_2` | Yes | R2 FASTQ; leave blank only for single-end Digenome data |
| `control` | No | Name of the matched control sample |
| `variant_vcf` | No | Indexed `.vcf.gz` containing known variants |

Important rules:

- controls are optional
- a named control must have its own row in the same samplesheet
- the control row leaves `control` blank
- one control may be shared by multiple treated samples
- control and VCF values must be consistent across rows for one sample
- each FASTQ path may appear only once in the samplesheet
- R1 and R2 cannot point to the same file
- unknown or misspelled columns are rejected
- FASTQs must end in `.fastq.gz` or `.fq.gz`
- nDigenome requires paired-end reads
- sample and control names use letters, digits, `.`, `_`, `-`

The template is at `assets/samplesheet_template.csv`.

## Configuration and parameters

For a complete lookup of all parameters, defaults, comparison semantics,
interactions, and examples, see [docs/parameters.md](docs/parameters.md).

Common runtime options:

| Option | Default | Purpose |
| --- | --- | --- |
| `--analysis` | `digenome` | Select DSB or SSB calling |
| `--genome` | Required | Configured genome name (e.g., `GRCh38`) |
| `--outdir` | `results` | Published results directory |
| `-work-dir` | Nextflow default | Temporary Nextflow work directory |
| `--keep_multimappers` | Off | Retain MAPQ 0 primary alignments for repetitive regions |
| `--cleavage_chunks` | `8` | Coordinate chunks per sample, maximum concurrent one-CPU tasks |
| `--genome_blacklist` | Unset | BED/BED.gz regions excluded before cleavage scanning |
| `--publish_trimmed_fastqs` | Off | Publish fastp output FASTQs |

Use `--cleavage_chunks 1` for serial cleavage calling. Increasing the value
allows more caller jobs to run concurrently but does not change the calling
rules or expected results. Large chromosomes are split into coordinate
intervals, so increasing the value can shorten a job that would otherwise
contain one whole chromosome. The planner may create fewer chunks only when
the mapped data and coordinate resolution cannot produce that many nonempty
ownership ranges.

Nextflow runtime options (one leading hyphen):

| Option | Default | Description |
| --- | --- | --- |
| `-profile` | None | Selects execution configuration: `standard`, `slurm`, `apptainer`, `longleaf` |
| `-resume` | Off | Reuses compatible cached tasks from the selected work directory |
| `-work-dir` | `work` | Nextflow task work directory; use stable, separate locations |
| `-c` | None | Adds a Nextflow configuration file (custom genomes, containers, resources) |

## Genomes

Configured genome names:

| Genome | Reference |
| --- | --- |
| `GRCh38` | `/proj/seq/data/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa` |
| `GRCh37` | `/proj/seq/data/HG19_UCSC/Sequence/WholeGenomeFasta/genome.fa` |
| `GRCm39` | `/proj/seq/data/GRCm39_GENCODE/Sequence/WholeGenomeFasta/GRCm39.primary_assembly.genome.fa` |

To use a different genome, add it in a custom configuration file:

```bash
nextflow run . -c my.config --genome MyGenome
```

where `my.config` contains:

```groovy
params {
    genomes {
        MyGenome = [fasta: '/path/to/genome.fa']
    }
}
```

## Genome blacklist

Use an optional genome-matched BED file to skip regions before endpoint
candidate scanning:

```bash
nextflow run . \
  ... \
  --genome_blacklist /path/to/excluded_regions.bed.gz
```

The file must use 0-based, half-open BED coordinates and BAM-matching contig
names. Overlapping or adjacent BED rows are merged. Unknown contigs,
out-of-range coordinates, malformed rows, and supplied files with no intervals
stop the run.

Focal nDigenome endpoints and both Digenome endpoints are excluded inside
blacklisted regions. Blacklisted opposite-strand endpoints also do not affect
nDigenome classification. The BED is optional and path-only: the pipeline does
not download, select, or maintain blacklist resources automatically.

The sample QC JSON records the staged filename, SHA-256 checksum, merged
interval count, and excluded-base count. Final filters and q-values remain
sample-wide over candidates from the callable, nonblacklisted genome.

## Execution profiles

| Profile | Use |
| --- | --- |
| `longleaf` | Longleaf with SLURM, Apptainer, site mounts, and shared paths |
| `slurm` | Generic SLURM execution |
| `apptainer` | Local execution with Apptainer |
| `standard` | Local configuration without a container engine |

Combine profiles with commas and the `test` profile last, e.g. `-profile apptainer,test` or
`-profile longleaf,test`.

## Multimappers

`--keep_multimappers` coordinates several settings:

- bwa-mem2 uses `-a` to retain secondary alignments
- MAPQ 0 primary alignments may contribute to calls
- the supporting-read mean-MAPQ filter is disabled
- fastp low-complexity filtering is disabled
- secondary and supplementary alignments remain diagnostic only
- duplicate-marked reads remain excluded

Each read counts once at BWA's selected primary location. It is not counted at
every reported secondary location. This prevents count inflation, but support
can be divided among equivalent repeat copies. The policy is recorded in
`analysis_parameters.json` and each sample QC JSON.

## Controls and variants

Controls are recommended but not required. Uncontrolled rows are labeled
`UNCONTROLLED` and are not filtered simply because a control is absent.

When a control is provided, the caller measures the same coordinate in the
control and calculates fold enrichment, Fisher exact p-values, and
Benjamini-Hochberg q-values. If local control depth is below
`--cleavage_control_min_depth` (default `1`), the row is labeled
`INSUFFICIENT_CONTROL_COVERAGE` and filtered without calculating misleading
fold, p, or q values.

The optional VCF must:

- be bgzip compressed as `.vcf.gz`
- have a `.tbi` or `.csi` index
- use the same reference build and contig names as the BAM
- contain every analyzed BAM contig

A mismatch such as `chr1` versus `1` stops the run rather than silently
disabling known-indel annotation.

## Calling, filtering, and output

Digenome mode pairs nearby forward and reverse 5-prime endpoint pileups.
nDigenome mode evaluates each strand independently and rejects meaningful
opposite-strand support from the high-confidence SSB set.

Digenome outputs report `digenome_pair_score`, which drives pairing and
filtering, plus `rgen_digenome_score`, which reproduces the standalone CRISPR
RGEN Tools v1.0 score for comparison with historical results. Matched controls
receive both scores independently.

Output tiers:

- `PASS`: retained in `*.high_confidence.tsv` and BED
- filtered rows with shared artifact evidence: retained in `*.artifact.tsv`
- other filtered rows: retained in `*.manual_review.tsv` for user review
- every row: retained in the complete `*.all.tsv` audit output

For the complete equations, exact comparisons, matching behavior, and all filter
reasons, see [docs/cleavage_algorithm.md](docs/cleavage_algorithm.md).

## Reference cache

bwa-mem2 indexes are stored under:

```text
<ref_cache>/<genome>/<fasta_sha256>/<bwa_mem2_version>/bwamem2/
```

By default the cache is relative to the checkout being run:

```text
<pipeline-directory>/reference_cache
```

The SHA-256 is calculated from the complete FASTA. A future run reuses an
index only when its genome name, FASTA hash, stable bwa-mem2 version, and
required index files match. Concurrent runs use a lock, so one run builds
while the other waits and then reuses the completed index.

The index helper sets shared-cache permissions explicitly instead of relying
on user umask or default ACL inheritance. Because the cache key includes both
the FASTA hash and bwa-mem2 version, each namespace directory through the
version level is mode `1777`:

```text
<ref_cache>/                         1777
<ref_cache>/<genome>/                1777
<ref_cache>/<genome>/<fasta_sha256>/ 1777
.../<bwa_mem2_version>/              1777
.../<bwa_mem2_version>/bwamem2/      755
.../bwamem2/index files              644
```

The sticky bit lets users create independent cache keys without deleting or
renaming entries owned by another user. Completed index directories are
globally readable/traversable but writable only by their creator. Build-lock
directories and owner records are also globally readable so concurrent users
can identify the active builder. Project filesystems may additionally preserve
setgid, producing modes such as `3777` or `2755`; the helper accepts that
additional group-inheritance bit.

## Outputs

Mode-specific outputs:

```text
<outdir>/<analysis>/
├── Sample.<analysis>.all.tsv
├── Sample.<analysis>.high_confidence.tsv
├── Sample.<analysis>.manual_review.tsv
├── Sample.<analysis>.artifact.tsv
├── Sample.<analysis>.bed
├── Sample.<analysis>.qc.json
└── Sample.<analysis>_mqc.tsv
```

Shared outputs:

```text
<outdir>/
├── bam/
├── fastp/
├── qc/
├── multiqc/
├── trimmed_fastqs/ (only with --publish_trimmed_fastqs)
└── pipeline_info/
```

`pipeline_info/analysis_parameters.json` records the resolved mode, thresholds,
reference, containers, multimapper policy, and chunk settings. The cleavage
chunk plan is stored under `pipeline_info/cleavage_chunks/`. It reports each
chunk's owned genomic intervals, callable and excluded bases, and estimated
total and callable mapped records.

## Testing

Run the local test suite:

```bash
python3 -m pip install -r requirements-test.txt
./tests/run_tests.sh
```

The suite uses Python `unittest` and synthetic indexed BAMs generated with
pysam. It covers workflow contracts, samplesheets, endpoint coordinates, DSB
pairing, SSB classification, controls, variants, artifacts, chunks, and
reference caching.

Before production use, run both modes through the Longleaf smoke test and
validate with known-positive material, untreated controls, and representative
library preparations. Repetitive-region work should include a relevant truth
set for `--keep_multimappers`.

See [docs/testing.md](docs/testing.md) for smoke-test commands and the limits
of local validation.
