# simgtsrch

C code to search for similar sets of genotypes.

E.g. we might have large number of accessions, and for each have
genotype information for a large number of markers.

Assumes an input file format like this:

MARKER mrkrid1 mrkrid2 ...
accid1  01101202011...
accid2  00101102011...
...

The genotypes are coded like this:
0 'reference' homozygous
1 heterozygous
2 alternative homozygous
any other (non-whitespace) character is treated as missing data.

All accessions must have a genotypes string of the same length.

Strategy is to compare genotype k-mers.
