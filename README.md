# simsearch and pedigree_test

///////////  simsearch:  //////////////////////

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



///////////  pedigree_test:  ////////////////////

C code to check pedigrees given in a file by comparing sets of genotypes.
In particular look at:

hgmr = (n02 + n21)/(n02 + n20 + n00 + n22)
this should be small when a parent-offspring relationship exists between the two genotypes.

and also at:
r = (n01 + n21)/(n01 + n21 + n00 + n22)

as a test of selfing. If, e.g. n01 is the number of genotypes pairs with the female parent's gt being 0, and the
offspring's gt being 1, then if in the self case, the male parent's gt will also be 0, and therefore the offspring's
gt should also be 0, not 1, so n01 (and n21) should be small, compared to n00 (n22), so r should be small. This only works if
we correctly know one parent, and are testing whether the other one is the same as that one, or distinct.

