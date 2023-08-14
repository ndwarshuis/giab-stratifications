# Stratification Validation

Prior to release, the stratifacations were validated according to the criteria
outlined here. The results of these tests are available in this directory for
reference.

## Unit tests

Each bed file is tested to meet the following criteria:
- must be a valid bed file with exactly 3 columns separated by tabs
- each chromosome must be present in the corresponding reference
- all lines must be sorted numerically by chromosome and start region (ie `X`
  and `Y` are 23 and 24 respectively)
- no regions can overlap
- all regions must be within chromosomal boundaries
- no regions can intersect with gaps as applicable (except for the gaps
  stratification itself obviously)
- each bed file is `bgzip`-ed

These tests do not produce any results so there is nothing to show here; the
pipeline that generated these stratifications runs these tests internally.

## Coverage plots

`coverage_plots.html` shows the coverage of each stratification file according
to the total length of each reference/chromosome. This was inspected prior to
release to ensure sensible coverage.

## Benchmarking

`benchmark_summary.html` contains benchmark results (using `hap.py`) for a
select subset of the stratifications. This was inspected prior to release for
sensible results.

## Comparison to previous versions

The directories `REF@all` contain comparisons between this version and the most
recent previous stratification version.

`diagnostics.tsv` contains a summary of the results of this comparison, and
contains columns as follows (note `A` means "this version" and `B` means
"previous version"):
- `total_A/B`: the total number of bases in the bed file
- `diff_AB`: the difference in total number of bases between `A` and `B`
- `total_shared`: the number of bases shared between `A` and `B`
- `shared_A/B`: the number of shared bases over the total number of bases in `A/B`
- `bedA/B`: the name of the bed file

Note that the names between `A` and `B` may differ slightly. Also, bed files
that are present in only `A` or `B` (not both) will have a blank in `bedA/B`.

The the `.bed.gz` files alongside `diagnostics.tsv` encode the exact differences
where they are present given a stratification level with both a current and
previous version. Each line represents a difference. These are regular bed files
with additional columns:
- `bed`: 0 or 1; 0 means the region in this line is in the current bed but not
  the previous bed (and vice versa)
- `adj`: an adjacency symbol indicating the relationship of the change in the file
  denoted by `bed` and the other bed:
  - `>`: the additional region is to the left of a region in the other bed
  - `<`: the additional region is to the right of a region in the other bed
  - `<>`: the additional region is in between two other regions in the other bed
  - `.`: the additional region does not border any regions in the other bed file
- `length`: end - start (for easy visualization)
- `other_bed`: the name of the previous bed (not important since this is in the
  filename)



