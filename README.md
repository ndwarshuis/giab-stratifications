# GIAB Stratifications Pipeline

This repository contains the snakemake pipeline for generating the
Genome-in-a-Bottle stratifications. This pipeline is designed to be reference
agnostic. For the complete stratifications generated for existing references
such as CHM13 or GRCh38, see
[here](https://github.com/ndwarshuis/giab-stratifications).

Note that this pipeline only applies to stratifications starting with version
v3.2. For older versions, please see
[here](https://github.com/genome-in-a-bottle/genome-stratifications) for the
equivalent methodological reference.

## Repository Overview

All stratifications can be generated using a snakemake pipeline from a reference
fasta file and accompanying set of source files.

### Included Categories

The following stratication categories can be generated programmatically by this
pipeline:
* Low Complexity
* Segmental Duplications
* XY
* Functional
* GC Content
* Mappability
* Telomeres
* Union

For other categories (such as ancestry) that may be derived from less
standardized means, it is possible to include raw bed files as-is. In practice,
this feature is often used for stratifications that are difficult to create and
that exist in previous versions (assuming they don't need to change).

### Requirements

#### Software

The only software requirement to install dependencies is `conda` (or `mamba`).

From the root of this repo, run the following:

```
mamba env create -f env.yml
```

This will create and env called `giab-strats`. Activate it with:

```
conda activate giab-strats
```

#### Hardware

If only generating stratifications for a few small chromosomes (as is the case
for testing and developing), the entire pipeline can be run in several minutes
on a consumer-grade laptop (8 cores, 16GB RAM). The only exception is generating
mappability (which runs GEM) which will require 10-20 minutes to complete for
even a small chromosome (like chr21).

TODO add cluster requirements

### Configuration

An example configuration can be found in `config` (currently used for testing
and continuous integration). This contains two files: `testing.yml` and
`testing-full.yml` which are used during integration testing and double as
example configurations. `testing.yml` is documented such so as to serve as a
starting point for users who wish to run this pipeline themselves.

Conceptually, the configuration is organized into two levels: `stratifications`
and `builds`. A `stratification` encodes the reference FASTA file along with
associated files used to generate the stratifications (which may be either
remote URLs are local files depending on what is available). Each
`stratification` has one or more `build`s which specify which stratifications
levels (eg LowComplexity, GCContent, etc) and chromosomes to include in the
final output.

### Running

To run the pipeline, use the following command (assuming the environment is set
up):

```
snakemake \
    --use-conda \
    -c <number_of_cores> --rerun-incomplete \
    --configfile=path/to/your/config/config.yml \
    all
```

Note that if any source files are specified as local in the configuration, they 
must exist or the pipeline will refuse to run.

### Output

All output will either be in `resources` (downloaded files) or `results`
(processed files). For most users the only important files will be in
`results/final`, which in turn has each build and reference specified as
`reference_key@build_key` (see configuration section and `config/testing.yml`
for more details on what these mean). Within each of these files are the
stratification BED files and associated metadata.
