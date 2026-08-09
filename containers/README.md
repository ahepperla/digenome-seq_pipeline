# Containers

The pipeline uses Nextflow's Apptainer integration. Process scripts do not
launch nested container runtimes.

## Image location

SIF files are excluded from Git because of their size. The shared Longleaf
installation stores the validated images in:

```text
/proj/jmsimon/Zylka/digenome-seq_pipeline/containers
```

Nextflow expects the files below in the `containers/` directory of whichever
pipeline checkout is being run. A separate checkout must have these images
copied or provisioned there before running the pipeline.

## Required images

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

`sources.tsv` records, for every SIF:

- the human-readable source tag
- the immutable multi-architecture OCI index digest
- the selected `linux/amd64` image-manifest digest
- whether the SIF was converted directly or built from a definition file
- the build architecture, Apptainer version, and UTC build time embedded in
  the SIF
- the UTC time when the tag and digest were last checked against Docker Hub

The source tags are readable aliases and may change upstream. The recorded OCI
digests are the immutable source identities. `checksums.sha256` separately
records the complete checksum of each built SIF.

The values in `sources.tsv` were checked against the metadata embedded in the
four SIFs and against Docker Hub on August 9, 2026. For publication or external
peer review, archive the exact SIF files in a persistent repository and cite
that accession or DOI. The shared Longleaf paths alone are not accessible to
outside reviewers.

The current cleavage SIF contains `pysam==0.23.3` and Debian `procps` version
`2:4.0.2-3`. Its definition pins the Python base-image digest and the pysam
version, but the original build did not pin a Debian repository snapshot or
the pysam wheel hash. The archived SIF checksum therefore identifies the exact
reviewed artifact; a later rebuild should be functionally equivalent but is
not guaranteed to have the same SIF checksum.
