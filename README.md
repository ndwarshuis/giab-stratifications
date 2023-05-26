# GIAB Stratifications

This the the pipeline that produces the official Genome-in-a-Bottle
stratifications, starting with version 3.2. All current and previous
stratifications can be found within the [NCBI FTP
site](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/).

Currently, this pipeline is set up to build stratifications for the following
references:
* GRCh37
* GRCh38
* CHM13

Further details (including all source files) can
be found in the configuration at `config/all.yml` in this repository as well
as the reference-agnostic [snakemake
pipeline](https://github.com/ndwarshuis/giab-strats-smk).

## Running this pipeline

This repository is primarily for documentation purposes, as the runtime is
defined in terms of NIST-specific resources.

However, for those wishing to run the pipeline themselves, here are several
options:

### No Cluster Environment

Just call the pipeline directly using snakemake:

```
snakemake --cores 20 --configfile=config/all.yml all
```

Note that most rules in the pipeline require 12MB of memory or less with the
exception of mappability (which runs GEM) and requires ~32GB for size of the
references here.

### Cluster Environment

First, create a profile from that located at
`workflow/profiles/nisaba/config.yml` (which is designed for a NIST-specific
cluster which runs Slurm).

Then run snakemake with this profile:

```
snakemake --cores 20 --configfile=config/all.yml --profile=workflow/profiles/yerprofile all
```
