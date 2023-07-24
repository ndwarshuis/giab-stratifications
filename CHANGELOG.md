# 3.0.0

- generalize chr prefix into a pattern (to allow recognizing chromosome names
  like "chr1_PATERNAL")
- relax constraints in input files; if not provided, output will not be
  generated; this is useful for cases where the input files do not exist.
  - affected strats: low complexity, xy, functional

# 2.10.0

- add AT/GC to low complexity just for homopolymers
- add AT/GC low complexity to benchmark output

# 2.9.0

- remove AT/GC from low complexity (for now)

# 2.8.1

- fix some random typos and bugs

# 2.8.0

- lower max low complexity length to 150bp

# 2.7.1

- fix typos

# 2.7.0

- make small rules not run on slurm
- fix lots of errors involving the final list of strats (actually involved using
  checkpoints for rules with complex output)
- fix formatting errors in checksum file (which didn't actually allow md5sun to
  run previously)

# 2.6.1

- fix chrom order in coverage plots

# 2.6.0

- add option to remove gaps from imported strat beds

# 2.5.1

- use plaintext and not gzipped file when performing md5 checks for resources

# 2.5.0

- make comparison way faster and more intuitive
- fix missing middle GC content file
- automatically make flanking gaps stratification

# 2.4.1

- fix overlapping level bug

# 2.4.0

- make benchmark subsets configurable

# 2.3.0

- automatically make gaps stratification

# 2.2.0

- add comparison functions to config/pipeline to test how much generated
  strats have changed relative to previous versions
- pipeline now fails on http 404 (or other bad request)

# 2.1.1

- fix typo

# 2.1.0

- make other strat name constraint more permissive

# 2.0.0

- add telomere stratification
- allow external beds to be used as stratifications
- add gc/at homopolymer bed files
- add benchmarking to validation postprocessing

# 1.0.0

- in the beginning there was darkness
