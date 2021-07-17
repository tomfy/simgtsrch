#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "vect.h"


typedef struct{
  Vchar* id;
  long index; // the index in the accessions array of GenotypesSet
  Vchar* genotypes;
  long missing_data_count;
}Accession;

typedef struct{
  long capacity; 
  long n_accessions; //
  long n_markers; //

  Accession* accessions; // array of accessions
  
  Vstr* accession_ids; // array of accession_ids
  Vlong* accession_missing_data_counts;
  Vstr* genotype_sets; //
  
  Vstr* marker_ids; // array of marker_ids
  Vlong* marker_missing_data_counts; //
  
}GenotypesSet;

Accession* construct_accession(char* id, long idx, char* genotypes);
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count);
void free_accession(Accession* the_accession);
void free_accession_innards(Accession* the_accession);

GenotypesSet* read_genotypes_file_and_store(FILE* g_stream, double delta, double max_missing_data_fraction);
void check_gtsset(GenotypesSet* gtsset);
GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts);
void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction);
GenotypesSet* construct_cleaned_genotypesset(GenotypesSet* the_gtsset, double max_md_fraction);
void set_accession_missing_data_counts(GenotypesSet* the_gtsset); 
void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset);
void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset);
void free_genotypesset(GenotypesSet* the_gtsset);
