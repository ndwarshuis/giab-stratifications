#remove all chroms that are NOT 1:22 X or Y.  Using removes chr prefix
gzcat hg38.gap.bed.gz | \
    sed 's/^chr//'| \
    grep -Ev '^M|^[0-9][0-9]_|[0-9]_|^Un' | \
    bgzip > hg38.gap_1thruY.bed

#change X and Y to number, sort, merge and sum.   Could add this to command above and convert 23/24 back to X and Y using sed 's/^[a-zA-Z0-9_]/chr&/'
gzcat hg38.gap_1thruY.bed.gz | \
    sed 's/^X/23/;s/^Y/24/'| \
     sort -k1,1n -k2,2n | \
     bedtools merge -i stdin -d 100| \
     awk '{ sum+=$3;sum-=$2 } END { print sum }'

# hg38.gap.bed.gz  219497085 1-22, X, Y + chrom contigs
# GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N.bed  165046090  1-22, X, Y + chrom contigs

# hg19.gap.bed.gz   239845127  1-22, X, Y + some with hap, e.g. chr6_apd_hap1, only for chroms 6 and 17
# example_of_no_ref_regions_input_file_b37.bed   200630157  1-22, X



# hg38.gap_1thruY.bed.gz   210157514   1-22, X, Y  includes centromeres
# GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N_1thuY.bed.gz  164556184  1-22, X, Y
# gap38_1thruY_sorted_merged.bed 150610728 1-22 X,Y

# hg19.gap_1thruY.bed.gz   234344806  1-22, X, Y
# hg19.gap_1thruX.bed.gz   200624806  1-22, X
# example_of_no_ref_regions_input_file_b37_1thruX.bed.gz   200630157  1-22, X


#Using trimmed (1:22, X and/or Y) create new N BEDs that are sorted and merged. 
gzcat hg38.gap_1thruY.bed.gz | \
     sed 's/^X/23/;s/^Y/24/' | \
     sort -k1,1n -k2,2n | \
     sed 's/^23/X/;s/^24/Y/' | \
    sed 's/^[a-zA-Z0-9_]/chr&/' | \
     bedtools merge -i stdin -d 100 > hg38.gap_1thruY_sorted_merged.bed
gzcat GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N_1thuY.bed.gz | \
     sed 's/^X/23/;s/^Y/24/'| \
     sort -k1,1n -k2,2n | \
     sed 's/^23/X/;s/^24/Y/' | \
     sed 's/^[a-zA-Z0-9_]/chr&/' | \
     bedtools merge -i stdin -d 100 > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N_1thruY_sorted_merged.bed

gzcat hg19.gap_1thruX.bed.gz | \
     sed 's/^X/23/;s/^Y/24/'| \
     sort -k1,1n -k2,2n | \
     sed 's/^23/X/;s/^24/Y/' | \
     bedtools merge -i stdin -d 100 > hg19.gap_1thruX_sorted_merged.bed
gzcat example_of_no_ref_regions_input_file_b37_1thruX.bed.gz | \
     sed 's/^X/23/;s/^Y/24/'| \
     sort -k1,1n -k2,2n | \
     sed 's/^23/X/;s/^24/Y/' | \
     bedtools merge -i stdin -d 100 > example_of_no_ref_regions_input_file_b37_1thruX_sorted_merged.bed
gzcat hg19.gap.bed.gz | \
     sed 's/^X/23/;s/^Y/24/'| \
     sort -k1,1n -k2,2n | \
     sed 's/^23/X/;s/^24/Y/' | \
     bedtools merge -i stdin -d 100 > hg19.gap_1thruY_sorted_merged.bed

#look at where the new N files differ
bedtools subtract -a hg38.gap_1thruY_sorted_merged.bed -b GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N_1thruY_sorted_merged.bed > subtractbed_38gapMINUS38refNs_1thruY.bed
bedtools subtract -a example_of_no_ref_regions_input_file_b37_1thruX_sorted_merged.bed  -b hg19.gap_1thruX_sorted_merged.bed > subtractbed_37refNsMINUS37gap_1thruX.bed
bedtools subtract -a GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N_1thruY_sorted_merged.bed -b gap_1thruY_sorted_merged.bed > subtractbed_38refNsMINUSnew38gap_1thruY.bed


What we did with N files -- all files in /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications
Noticed discrepency in validation where BED coverage > 1.  The thought was that this was 
because Ns were included in stratification beds however we had removed them from reference.
We decided to remove Ns from stratification beds and again look at coverage using 
-GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_REF_N.bed
-example_of_no_ref_regions_input_file_b37.bed.  
This improved coverage however Y for ref37 was not present in example_of_no_ref_regions_input_file_b37.bed
We then decided to try to use Aarons N beds
-hg38.gap.bed.gz
-hg19.gap.bed.gz
I first trimmed all (Aaron and one JW provided) files to just include the chroms we were 
interested in (1:22, XY). I then sorted and merged these files.  We looked at total region
size covered by each and compared files as appropriate.  I perfomed subract bed so JZ
could look specifically at where the regions differed.  37Ns for both were pretty much
the same with exception of 5K bases.  Aarons gap 38 file was significantly larger than
the one JW provided an upon inspection JZ found many centeromere regions were covered by 
aarons bed.  The README for Aarons gap bed appears to show addition of centromeres to gap
file (see below)
## Gap+centromere track

LTMPDIR=$(mktemp -d)
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz $LTMPDIR
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz $LTMPDIR
zcat $LTMPDIR/gap.txt.gz | \
    cut -f2-4,8 | \
    awk '{n=$4; if(n=="clone" || n=="contig" || n=="scaffold" || n=="short_arm") { n="gap"; } print $1"\t"$2"\t"$3"\t"n; }' | \
    cat - <(zcat $LTMPDIR/centromeres.txt.gz | awk '{ print $2"\t"$3"\t"$4"\tcentromere"; }') | \
    sort -k1,1 -k2,2g | bgzip -c > annotation/gap.bed.gz
tabix -f annotation/gap.bed.gz
rm -rf $LTMPDIR && unset LTMPDIR

## Centromere track
LTMPDIR=$(mktemp -d)
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz $LTMPDIR
zcat $LTMPDIR/gap.txt.gz | \
    cut -f2-4,8 | \
    awk '{n=$4; if(n=="clone" || n=="contig" || n=="scaffold" || n=="short_arm") { n="gap"; } print $1"\t"$2"\t"$3"\t"n; }' | \
    sort -k1,1 -k2,2g | \
    bgzip -c > annotation/gap.bed.gz
tabix -f annotation/gap.bed.gz
rm -rf $LTMPDIR && unset LTMPDIR

The "centromere track" appears to be what we want to I ran the code
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz .
gzcat gap.txt.gz | \
    cut -f2-4,8 | \
    awk '{n=$4; if(n=="clone" || n=="contig" || n=="scaffold" || n=="short_arm") { n="gap"; } print $1"\t"$2"\t"$3"\t"n; }' |
    sort -k1,1 -k2,2g | \
    bgzip -c > gap38.bed.gz
tabix -f gap38.bed.gz

looking at this output there are no "centromeres" however there are in the orginal gap
file we have for Aaron (hg38.gap.bed.gz)

gzcat gap38.bed.gz | \
    sed 's/^chr//' | \
    grep -Ev '^M|^[0-9][0-9]_|[0-9]_|^Un' | \
    sed 's/^X/23/;s/^Y/24/' | \
    sort -k1,1n -k2,2n | \
    sed 's/^23/X/;s/^24/Y/' | \
    sed 's/^[a-zA-Z0-9_]/chr&/' | \
    bedtools merge -i stdin -d 100 > gap_38_1thruY_sorted_merged.bed

I ran aarons script for 37 gap just to confirm there were no differences with his file and
they matched.

Moving forward for Ns we will use:
37 --> hg19.gap_1thruY.bed.gz
38 --> gap38_1thruY_sorted_merged.bed

# 12/10/19
# Created psuedoautosomal region bed files (PSA_Y_GRCh38.bed
# PSA_Y_hg19.bed) for Y using information found on UCSC genome site 
# PSA regions for Y (GRCh38) : chrY:10,000-2,781,479 and chrY:56,887,902-57,217,415
# PSA regions for Y (hg19) : chrY:10,001-2,649,520 and chrY:59,034,050-59,363,566
# These regions really should be Ns so we are going to subtract them from the genome BED files


# JZ and Nate want to use the same methods to subtract out Ns for the genomes as we did 
# for the stratifications.  In lieu of Nate's R function to create beds we will do the 
# following
# 1) Use the following genome beds from JZ (slack message 12/10/19)
# human.b37.1_22XY.genome.bed
# human.hg38.1_22XY.genome.bed

sort the genome BEDs
sed 's/^X/23/;s/^Y/24/' /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.b37.1_22XY.genome.bed | \
sort -k1,1n -k2,2n  | \
sed 's/^23/X/;s/^24/Y/' > /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.b37.1_22XY.genome.sorted.bed

sed 's/^chr//' /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.hg38.1_22XY.genome.bed | \
sed 's/^X/23/;s/^Y/24/' |\
sort -k1,1n -k2,2n  | \
sed 's/^23/X/;s/^24/Y/' | \
sed 's/^[a-zA-Z0-9_]/chr&/' > /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.hg38.1_22XY.genome.sorted.bed

# 2) Remove N's in reference. Subtract gap BEDs from genome BEDs
subtractBed -a /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.b37.1_22XY.genome.sorted.bed -b /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/refNs/hg19.gap_1thruY_sorted_merged.bed > /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.b37.1_22XY.genome.sorted_noNs.bed
subtractBed -a /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.hg38.1_22XY.genome.sorted.bed -b /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/refNs/gap_38_1thruY_sorted_merged.bed > /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.hg38.1_22XY.genome.sorted_noNs.bed

# 3) Remove psuedoautosomal regions for Y in reference BED using subtractBed 
subtractBed -a /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.b37.1_22XY.genome.sorted_noNs.bed -b /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/refNs/PSA_Y_hg19.bed > /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.b37.1_22XY.genome.sorted_noNs_noPSAY.bed
subtractBed -a /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.hg38.1_22XY.genome.sorted_noNs.bed -b /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/refNs/PSA_Y_GRCh38.bed > /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/reference_beds/human.hg38.1_22XY.genome.sorted_noNs_noPSAY.bed

# 4) Merge and sum by chrom in genome BED use this as denominator for calculating coverage
mergeBed -i human.b37.1_22XY.genome.sorted_noNs_noPSAY.bed -d 100 > human.b37.1_22XY.genome.sorted_noNs_noPSAY_merged.bed
mergeBed -i human.hg38.1_22XY.genome.sorted_noNs_noPSAY.bed -d 100 > human.hg38.1_22XY.genome.sorted_noNs_noPSAY_merged.bed
