# Test files for CHM13 Low Complexity

The following are reduced files (to keep repo size down) for testing
LowComplexity.

## TRF

Source file was downloaded from
[here](https://app.globus.org/file-manager?origin_id=9db1f0a6-a05a-11ea-8f06-0a21f750d19b&origin_path=%2Fteam-segdups%2FAssembly_analysis%2FMasked%2F/T2T_CHM13v2_trf.bed).

Source was filtered for chr21 and 22 and all but first three columns were
stripped:

```
gunzip -c src_trf.bed.gz | \
    grep -E '^(chr21|chr22)' | \
    cut -f1-3 | \
    gzip -c > trf21and22.bed.gz
```

## Censat

Source file from [slack
thread](https://t2t-consortium.slack.com/files/ULT7E06GL/F039A96RY84/t2t_censat_chm13v2.0_trackv2.0.bed).

Source was filtered for 21 and 22 (and the header was kept to keep it more
realistic):

```
gunzip -c src_censat.bed.gz | \
    grep -E '^(track|chr21|chr22)' | \
    gzip -c > censat21and22.bed.gz
```

