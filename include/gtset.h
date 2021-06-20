#include <stdlib.h>
#include <stdio.h>
#include "vect.h"


typedef struct{
  long n_accessions; //
  long n_markers; //
  Vstr* accession_ids; // array of accession_ids
  Vstr* genotype_sets; //
  Vstr* marker_ids; // array of marker_ids
  Vlong* marker_missing_data_counts; //
}GenotypesSet;

GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts);
void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction);
GenotypesSet* construct_cleaned_genotypesset(GenotypesSet* the_gtsset, double max_md_fraction);
void print_genotypesset(GenotypesSet* the_gtsset);
void print_genotypesset_summary_info(GenotypesSet* the_gtsset);
void free_genotypesset(GenotypesSet* the_gtsset);
