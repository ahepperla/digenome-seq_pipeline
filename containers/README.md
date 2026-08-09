# Containers

The pipeline uses Nextflow's Apptainer integration. Process scripts do not
launch nested container runtimes.

## Bundled images

| File | Purpose |
| --- | --- |
| `fastp_v1.3.3.sif` | read trimming |
| `bwa-mem2_v2.3_samtools_v1.22.sif` | alignment and BAM processing |
| `multiqc_v1.35.sif` | samplesheet validation and MultiQC |
| `cleavage_pysam_v0.23.3.sif` | unified Digenome/nDigenome caller |

## Cleavage image

`cleavage_pysam_v0.23.3.def` pins:

- Python 3.12.11 slim-bookworm by immutable registry digest
- `pysam==0.23.3`
- `procps`, which provides `ps` for Nextflow task metrics

Rebuild this image on a host with Apptainer:

```bash
cd /path/to/digenome-seq_pipeline
./containers/build_cleavage.sh
```

The build writes:

```text
cleavage_pysam_v0.23.3.sif
```

and updates its entry in `checksums.sha256`.
It also verifies the pinned pysam version and the presence of `ps`.

## Verification

Run these commands from the `containers` directory:

```bash
apptainer exec fastp_v1.3.3.sif fastp --version
apptainer exec bwa-mem2_v2.3_samtools_v1.22.sif bwa-mem2 version
apptainer exec bwa-mem2_v2.3_samtools_v1.22.sif samtools --version
apptainer exec multiqc_v1.35.sif multiqc --version
apptainer exec cleavage_pysam_v0.23.3.sif \
  python3 -c 'import importlib.metadata; print(importlib.metadata.version("pysam"))'
sha256sum -c checksums.sha256
```

## Provenance

`sources.tsv` records immutable base-image digests. `checksums.sha256` records
complete SIF checksums.
