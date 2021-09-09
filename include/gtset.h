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
  long size;
  Accession** a;
}Vaccession;

typedef struct{
  long capacity;
  double delta;
  double max_marker_missing_data_fraction;
  long n_accessions; //
  long n_markers; //

  Vaccession* accessions;
  
  Vstr* marker_ids; // array of marker_ids
  Vlong* marker_missing_data_counts; //
  
}GenotypesSet;

// *****  functions  *****
// *****  Accession  *****
Accession* construct_accession(char* id, long idx, char* genotypes, long accession_md_count);
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count);
void free_accession(Accession* the_accession);
void free_accession_innards(Accession* the_accession);

// *****  Vaccession  *****
Vaccession* construct_vaccession(long cap);
void add_accession_to_vaccession(Vaccession* the_vacc, Accession* the_acc);
void free_vaccession(Vaccession* the_vacc);

// *****  GenotypesSet  *****
GenotypesSet* read_dosages_file_and_store(FILE* g_stream, double delta);
GenotypesSet* read_genotypes_file_and_store(FILE* g_stream);
void check_gtsset(GenotypesSet* gtsset);
GenotypesSet* construct_genotypesset(Vaccession* accessions, Vstr* marker_ids, Vlong* md_counts, double delta, double max_marker_md_fraction);
void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction);
GenotypesSet* construct_cleaned_genotypesset(const GenotypesSet* the_gtsset, double max_md_fraction);
void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset);
void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset);
void free_genotypesset(GenotypesSet* the_gtsset);

Vidxid* construct_vidxid(const GenotypesSet* the_gtsset);
Vidxid* construct_sorted_vidxid(const GenotypesSet* the_gtsset);
long check_idxid_map(Vidxid* vidxid, const GenotypesSet* the_gtsset);
