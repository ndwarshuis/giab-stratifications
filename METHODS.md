# Overview

The following is an overview of the pipeline, including data sources and
computational methods for each stratification type.

NOTE: everything in this document applies to a fully specified configuration (ie
all options in the `include` block are set to `true` and all chromosomes are
included). If this is not the case, some of these will not be applicable.

# 1. Stratification Generation

The following is a description of the generation process of each stratification
type up until postprocessing (see next).

## Functional Regions

Coding regions are based on the GFF and FTBL files from the RefSeq database. The
pipeline filters the GFF dataframe for type `CDS` and joins this to the FTBL
dataframe using the accession number. Chromosome coordinates are selected from
this to make a bed file.

## GC Content

The [seqtk](https://github.com/lh3/seqtk) tool is used to identify up to 100bp
regions with greater than/less than some percentage of GC content. The pipeline
generates multiple bed files at at these known percentage cutoffs, using the
fasta reference as the only input file. After adding 50bp slop to each region
and merging, adjacent cutoffs are subtracted from each other to produce banded
region files (for example, subtracting GC <15% from GC <20% gives all regions
between 15-20% GC). These banded files are part of the final stratifiation bed
file set of this type.

In addition, the pipeline generates beds for GC <25% / >65% and GC <30% / >55%.

The decision to include each of these bands was determined from internal
exploratory work as well as past data regarding high error rates due to
sequencing biases ([ref](doi:10.1186/gb-2013-14-5-r51)).

## Mappability

Hard to map regions are found using
[GEM](https://doi.org/10.1371/journal.pone.0030377), which is a program that
aligns sequences to themselves to find regions that do not align anywhere else
given a set of alignment parameters. These parameters include `L` (the length of
the regions to compare), `M` (the number of mismatches, ie SNVs), and `E` (the
number of gaps, ie INDELs).

The parameters that will be included in these stratifications are fully
configurable; however, a typical set is in the ballpark of `L`=100, `M`=2, and
`E`=1 (low stringency) or `L`=250, with `M`/`E` = 0 (high stringency).

NOTE: this is the most computationally heavy part of the entire pipeline. For
reference, it takes about an 30 minutes to run the indexer and 60 minutes to run
mappability with the high stringency parameters, both with 32 cores and 12 GB
memory. By way of comparison, most rules in the pipeline can run in several
minutes or less with low/constant memory (since they are just piped bedtools
commands) using a single thread.

## Low Complexity

In general, this type of stratification type combines different flavors of
repeats:

* **Simple Repeats**: tandem repeats as found from the TRF tool (for common
  assemblies available in the UCSC Genome Browser)
* **Repeat Masker**: repeats as defined using the rmsk tool (again usually
  available through UCSC genome browser). This pipeline only uses the
  `Low_complexity`, `Simple_repeat`, and `Satellite` classes (although see below
  for exception).
* **Censat**: if provided, the pipeline will use an alternative satellite bed
  files instead of the `Satellite` class from repeat masker (relevant for T2T).
* **Uniform Repeats**: unlike the above, these repeats are generated on-the-fly
  using a program (described further below).
  
### Uniform Repeats

These come in two subflavors: perfect and imperfect.

A "perfect uniform repeat" is defined as a repeating pattern of either 1, 2, 3,
or 4bp motifs (denoted lexically as "homopolymer", "di", "tri", and "quad"
repeats) that are less than or equal to a given length. These are generated
using a [C program](https://github.com/ndwarshuis/repseq) which scans the
reference for such repeats, producing a bed file annotated with the motif in the
4th column.

The "imperfect" analogue is any number of "perfect uniform repeats" separated by
1 non-repeat base (actually produced using `mergeBed -d 1` on any of the perfect
uniform repeat beds). Coloquially, a special case of these where the repeat
motif is 1 are known as "imperfect homopolymers" (eg long homopolymer stretches
with a few non-homopolymer bases scattered throughout).

## Segmental Duplications

Segmental duplications are created from the GenomicSuperDups database provided
by UCSC. Regions 100 bp apart or less are merged. The pipeline will also produce
separate files for segdups >10000bp long.

## XY

These stratification files are generated failry straightforwardly with minimal
modification to source data. PAR1/2 regions are specied directly in the
configuration file for both X and Y and used as is. Ampliconic/XTR regions are
used drectly from their source bed files (which can be obtained from [Melissa
Wilson's](https://github.com/SexChrLab/SexChrCoordinates) lab for common
genomes).

## Union

These are sevearl overall combinations of bedfiles from other types.

Currently the pipeline will combine segmental duplications and mappability
stratifications (eg all "hard to map" regions), and will further combine the
previous two types with GC content <30% / >65%, all low complexity regions, and
all XY features (XTR + Ampliconic depending on what is available).

# 2. Post Processing

## Merging and gap removal

All bed files are merged using `mergeBed` (with `-d` set to 0) to simplify all
overlapping regions. If a gaps bed file is specified, these gaps are removed
from each stratification bed file using `intersectBed`.

## Unit Testing

Each final stratification bed file is unit tested to ensure the following
properties:

* each file is in `bgzip` format
* each bed file has 3 columns, the first of which is a valid chromosome
  identifier and the second/third of which are positive integers
* the bed file is sorted
* the bed file has no overlapping regions
* the regions in the bed file are in valid range (not in gaps and not beyond
  chromosome boundaries)
  
## Qualitative Validation

Each stratification build comes with an autogenerated html file depicting the
coverage of each stratification bed file for each chromosome. These are to be
manually inspected as a sanity check, and are available to end users for their
scrutiny.
