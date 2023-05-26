# GIAB Stratifications

## Introduction

This repository contains the files used for generating and evaluating the v3.1
stratifications. The stratification files were developed with the Global
Alliance for Genomic Health (GA4GH) Benchmarking Team, the [Genome in a Bottle
Consortium (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) and the
[Telomere-to-Telomere Consortium
(T2T)](https://sites.google.com/ucsc.edu/t2tworkinggroup). They are intended as
a standard resource of BED files for use in stratifying true positive, false
positive and false negative variant calls in challenging and targeted regions of
the the genome.

These files can be used as a standard resource of BED files for use with GA4GH
benchmark tools such as [hap.py](https://github.com/Illumina/hap.py).

*NOTE: stratification BED files are only accessible on the [GIAB FTP
site](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/)*
 
## Summary Of Stratifications

Stratifications can be binned into several types: Low Complexity, Functional
Technically Difficult, Genome Specific, Functional Regions, GC content,
mappability, Other Difficult, Segmental Duplications, Union, Ancestry and XY.
General information for stratification types are provided below.

### Low Complexity

Regions with different types and sizes of low complexity sequence, e.g.,
homopolymers, STRs, VNTRs and other locally repetitive sequences.

### Segmental Duplications

Regions with segmental duplications (generally defined as repeated regions >1kb
with >90% similarity).

### XY

Chomosome XY specific regions such as PAR, XTR or ampliconic.

### Functional Regions

Regions to stratify variants inside and outside coding regions.

### GC Content

Regions with different ranges (%) of GC content.

### Mappability

Regions where short read mapping can be challenging.

### Telomeres

Telomeric regions

### Union

Regions with different general types of difficult regions or any type of
difficult region or complex variant. For example, performance can be measured in
just "easy" or "all difficult" regions of the genome.

### Not yet implemented

These stratifications existed in v3.1 but have yet to be implemented in
the current version. If these are desired, the configuration allows including
pre-constructed bed files as stratifications, so these can be obtained from
previous versions on new builds.

#### Ancestry

Regions with inferred patterns of local ancestry.

#### Genome Specific (GIAB benchmark v4.2.1)

Difficult regions due to potentially difficult variation in a NIST/GIAB sample,
including 1) regions containing putative compound heterozygous variants 2) small
regions containing multiple phased variants, 3) regions with potential
structural or copy number variation.

#### Functional Technically Difficult

Functional, or potentially functional, regions that are also likely to be
technically difficult to sequences.

#### Other Difficult

Highly variable regions like the VDJ and MHC, near gaps in the reference or
errors in the reference and rDNA (CHM13 only).

## Repository Overview

All stratifications can be generated using a snakemake pipeline from a reference
fasta file and accompanying set of source files.

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

## General Information

Author Information
- Principal Investigator: Justin Zook, NIST, jzook@nist.gov
- Nate Olson, NIST, nathanael.olson@nist.gov
- Justin Wagner, NIST, justin.wagner@nist.gov
- Jennifer McDaniel, NIST, jennifer.mcdaniel@nist.gov
- Nate Dwarshuis, NIST, nathan.dwarshuis@nist.gov

## Sharing/Access Information

Licenses/restrictions placed on the data, or limitations of reuse:
Publicly released data are freely available for reuse without embargo.

Citations for stratifications are located in the associated READMEs.

If stratifications were used in benchmarking with GA4GH/GIAB best practices or
hap.py please reference:

	Krusche, P., Trigg, L., Boutros, P.C. et al. 
	Best practices for benchmarking germline small-variant calls in human genomes. 
	Nat Biotechnol 37, 555-560 (2019). https://doi.org/10.1038/s41587-019-0054-x

### Links to publicly accessible locations of the data:

[GIAB FTP URL](https://ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/)
- Individual stratification BED files as well as zipped directories (tar.gz) of files
- stratification READMEs
- .tsvs for benchmarking with hap.py
- MD5 checksums

## Data-Use Policy 

This data/work was created by employees of the National Institute of Standards
and Technology (NIST), an agency of the Federal Government. Pursuant to title 17
United States Code Section 105, works of NIST employees are not subject to
copyright protection in the United States. This data/work may be subject to
foreign copyright.

The data/work is provided by NIST as a public service and is expressly provided
AS IS. NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR
STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA
ACCURACY. NIST does not warrant or make any representations regarding the use of
the data or the results thereof, including but not limited to the correctness,
accuracy, reliability or usefulness of the data. NIST SHALL NOT BE LIABLE AND
YOU HEREBY RELEASE NIST FROM LIABILITY FOR ANY INDIRECT, CONSEQUENTIAL, SPECIAL,
OR INCIDENTAL DAMAGES (INCLUDING DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS
INTERRUPTION, LOSS OF BUSINESS INFORMATION, AND THE LIKE), WHETHER ARISING IN
TORT, CONTRACT, OR OTHERWISE, ARISING FROM OR RELATING TO THE DATA (OR THE USE
OF OR INABILITY TO USE THIS DATA), EVEN IF NIST HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

To the extent that NIST may hold copyright in countries other than the United
States, you are hereby granted the non-exclusive irrevocable and unconditional
right to print, publish, prepare derivative works and distribute the NIST data,
in any medium, or authorize others to do so on your behalf, on a royalty-free
basis throughout the world.

You may improve, modify, and create derivative works of the data or any portion
of the data, and you may copy and distribute such modifications or works.
Modified works should carry a notice stating that you changed the data and
should note the date and nature of any such change. Please explicitly
acknowledge the National Institute of Standards and Technology as the source of
the data: Data citation recommendations are provided at
https://www.nist.gov/open/license.

Permission to use this data is contingent upon your acceptance of the terms of
this agreement and upon your providing appropriate acknowledgments of NISTâ€™s
creation of the data/work.
