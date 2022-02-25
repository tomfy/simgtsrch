#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "vect.h"

#define DO_ASSERT 1
#define UNKNOWN -1
#define DOSAGES 0
#define GENOTYPES 1
#define MISSING_DATA_CHAR '-'  
#define MISSING_DATA_VALUE -1
#define INIT_VACC_CAPACITY 2000

typedef struct{
  Vchar* id;
  long index; // the index in the accessions array of GenotypesSet
  // long n_markers;
  Vchar* genotypes;
  Vlong* chunk_patterns;
  long md_chunk_count;
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
  long n_accessions; // redundant.
  long n_ref_accessions; 
  long n_markers; // redundant.

  Vaccession* accessions;
  
  Vstr* marker_ids; // array of marker_ids
  Vlong* marker_missing_data_counts; //
  
}GenotypesSet;

// *****  functions  *****
long int_power(long base, long power);
long determine_file_format(char* filename);

// *****  Accession  *****
Accession* construct_accession(char* id, long idx, char* genotypes, long accession_md_count);
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count);

Accession* construct_accession_a(char* id, char* gtset);
long set_accession_chunk_patterns(Accession* the_gts, Vlong* m_indices, long n_chunks, long k);
char* print_accession(Accession* the_gts, FILE* ostream);

void free_accession(Accession* the_accession);
void free_accession_innards(Accession* the_accession);

// *****  Vaccession  *****
Vaccession* construct_vaccession(long cap);
void add_accession_to_vaccession(Vaccession* the_vacc, Accession* the_acc);
void set_vaccession_chunk_patterns(Vaccession* the_accessions, Vlong* m_indices, long n_chunks, long k);
void print_vaccession(Vaccession* the_accessions, FILE* ostream);
void check_accession_indices(Vaccession* the_accessions);
void free_vaccession(Vaccession* the_vacc);

// *****  GenotypesSet  *****
GenotypesSet* construct_empty_genotypesset(double delta, double max_marker_md_fraction);
GenotypesSet* read_dosages_file_and_store(char* input_filename, double delta);
GenotypesSet* read_genotypes_file_and_store(char* input_filename);

void add_accessions_to_genotypesset_from_file(char* input_filename, GenotypesSet* the_genotypes_set);
void read_dosages_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set);
void read_genotypes_file_and_add_to_genotypesset(char* input_filename, GenotypesSet* the_genotypes_set);

void check_gtsset(GenotypesSet* gtsset);
GenotypesSet* construct_genotypesset(Vaccession* accessions, Vstr* marker_ids, Vlong* md_counts, double delta, double max_marker_md_fraction);
void check_genotypesset(GenotypesSet* gtss);
GenotypesSet* construct_cleaned_genotypesset(const GenotypesSet* the_gtsset, double max_md_fraction);
void clean_genotypesset(GenotypesSet* the_genotypes_set);
void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset);
void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset);
void free_genotypesset(GenotypesSet* the_gtsset);

Vidxid* construct_vidxid(const GenotypesSet* the_gtsset);
Vidxid* construct_sorted_vidxid(const GenotypesSet* the_gtsset);
long check_idxid_map(Vidxid* vidxid, const GenotypesSet* the_gtsset);
