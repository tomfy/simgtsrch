// C version of program to test pedigrees using genotype data.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include "vect.h"

#define DBUG 1

  
typedef struct{
  long n_accessions; //
  long n_markers; //
  Vstr* accession_ids; // array of accession_ids
  Vstr* genotype_sets; //
  Vstr* marker_ids; // array of marker_ids
  Vlong* marker_missing_data_counts; //
}GenotypesSet;

GenotypesSet* read_genotypes_file_and_store(FILE* g_stream, double delta, double max_missing_data_fraction);

GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  if(gsets->size != acc_ids->size){ fprintf(stderr, "Inconsistency in construct_genotypesset\n"); exit(EXIT_FAILURE); }
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  the_gtsset->n_accessions = acc_ids->size; 
  the_gtsset->n_markers = marker_ids->size;
  the_gtsset->accession_ids = acc_ids;
  the_gtsset->marker_ids = marker_ids;
  the_gtsset->genotype_sets = gsets;
  the_gtsset->marker_missing_data_counts = md_counts;
  return the_gtsset;
}

void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction){
  assert(gtss->accession_ids->size == gtss->n_accessions);
  assert(gtss->genotype_sets->size == gtss->n_accessions);
  assert(gtss->marker_ids->size == gtss->n_markers);
  long* md_counts = (long*)calloc(gtss->n_markers, sizeof(long)); 
  for(long i=0; i<gtss->genotype_sets->size; i++){
    assert(strlen(gtss->genotype_sets->a[i]) == gtss->n_markers);
    for(long j=0; j<gtss->n_markers; j++){
      if(gtss->genotype_sets->a[i][j] == '3')md_counts[j]++;
    }
  }
  for(long j=0; j<gtss->n_markers; j++){
    assert(md_counts[j] == gtss->marker_missing_data_counts->a[j]);
  }
  free(md_counts);
  fprintf(stderr, "Successfully completed check_genotypesset\n");
}

GenotypesSet* construct_cleaned_genotypesset(GenotypesSet* the_gtsset, double max_md_fraction){
 
  Vstr* gsets = the_gtsset->genotype_sets;
  Vlong* md_counts = the_gtsset->marker_missing_data_counts;

  // identify the markers to keep:
  long n_markers_to_keep = 0;
  Vlong* md_ok = construct_vlong_zeroes(md_counts->size);
  Vstr* cleaned_marker_ids = construct_vstr(1000);
  Vlong* cleaned_md_counts = construct_vlong(1000);
  for(long i=0; i<md_counts->size; i++){
    if(md_counts->a[i] <= max_md_fraction*the_gtsset->n_markers){
      md_ok->a[i] = 1;
      n_markers_to_keep++;
      long marker_id_length = strlen(the_gtsset->marker_ids->a[i]);
      char* marker_id_to_keep = strcpy((char*)malloc((marker_id_length+1)*sizeof(char)), the_gtsset->marker_ids->a[i]);
      add_string_to_vstr(cleaned_marker_ids, marker_id_to_keep);
      add_long_to_vlong(cleaned_md_counts, md_counts->a[i]);
    }
  }
  fprintf(stderr, "after 1st loop\n");
  //GenotypesSet* the_cleaned_gtsset = (GenotypesSet*)malloc(sizeof(GenotypesSet));
  Vstr* copy_of_accids = construct_vstr_copy(the_gtsset->accession_ids);
  Vstr* cleaned_gsets = construct_vstr(the_gtsset->n_accessions);
  for(long i=0; i<gsets->size; i++){ // loop over accessions
    char* raw_gts = gsets->a[i]; // the string with all the genotypes for accession i
    char* cleaned_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    // char* cleaned_marker_ids = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      if(md_ok->a[j] == 1){
	cleaned_gts[k] = raw_gts[j];
	//	long len = strlen(the_gtsset->marker_ids[j]);
	// cleaned_marker_ids[k] = strcpy((char*)malloc((len+1)*sizeof(char), the_gtsset->marker_ids[j]);
	k++;
      }
    }
    cleaned_gts[k] = '\0'; // terminate with null.
    add_string_to_vstr(cleaned_gsets, cleaned_gts);
    if(DBUG) assert(k == n_markers_to_keep);
  }
  free_vlong(md_ok);
  //  GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  GenotypesSet* cleaned_gtsset = construct_genotypesset(copy_of_accids, cleaned_marker_ids, cleaned_gsets, cleaned_md_counts);
  return cleaned_gtsset;
}


void print_genotypesset(GenotypesSet* the_gtsset){
  printf("MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    printf("%s ", the_gtsset->marker_ids->a[i]);
  }printf("\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    printf("%s  %s\n", the_gtsset->accession_ids->a[i], the_gtsset->genotype_sets->a[i]);
  }
}

void print_genotypesset_summary_info(GenotypesSet* the_gtsset){
  fprintf(stderr, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  fprintf(stderr, "# n_markers: %ld\n", the_gtsset->n_markers);
}

void free_genotypesset(GenotypesSet* the_gtsset){
  free_vstr(the_gtsset->accession_ids);
   free_vstr(the_gtsset->marker_ids);
  free_vstr(the_gtsset->genotype_sets);
  free_vlong(the_gtsset->marker_missing_data_counts);
  free(the_gtsset);
}

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  double delta = 0.05; // default; control this with -d command line option.
  double max_marker_missing_data = 0.2; // default; control this with -x command line option.
  long max_number_of_accessions = 1000000;

  
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage:  %s -g <genotypes_file>  -p <pedigree_file>  options -d -x \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* genotypes_filename = NULL;
  FILE *g_stream = NULL;
  char* pedigrees_filename = NULL;
  FILE *p_stream = NULL;

  // g: genotypes filename, p pedigree filename,  d (delta for rounding), x max fraction of missing data for markers.
  int c;
  while((c = getopt(argc, argv, "g:p:d:x:")) != -1){
    // fprintf(stderr, "c: %c %s %d\n", c, optarg, optind);
    switch(c){
    case 'g':
      genotypes_filename = optarg;
      g_stream = fopen(genotypes_filename, "r");
      if(g_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	exit(EXIT_FAILURE);
      }
      break;
    case 'p':
      pedigrees_filename = optarg;
      p_stream = fopen(pedigrees_filename, "r");
      if(p_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", pedigrees_filename);
	exit(EXIT_FAILURE);
      }
      break;
    case 'd':
      if(optarg == 0){
	perror("option d requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	delta = atof(optarg);
      }
      break;
    case 'x':
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_marker_missing_data = atof(optarg);
      }
      break;
    case '?':
      printf("? case in command line processing switch.\n");
      if ((optopt == 'i') || (optopt == 'n') || (optopt == 'k'))
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf(stderr, "Unknown option character: %d\n", optopt);
      exit(EXIT_FAILURE);
    default:
      perror("default case (abort)\n");
      abort ();
    } // end of switch block
  } // end of loop over c.l. arguments
    // printf("optind: %d argc: %d\n", optind, argc);
  if(optind < argc){
    perror("Non-option arguments. Bye.\n");
    exit(EXIT_FAILURE);
  }
  if(genotypes_filename == NULL){
    perror("must specify genotype filename: -i <filename>");
    exit(EXIT_FAILURE);
  }
  if(pedigrees_filename == NULL){
    perror("must specify pedigrees filename: -i <filename>");
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "# genotypes file: %s  pedigree file: %s  delta: %5.3lf  max marker missing data: %5.3lf \n",
	  genotypes_filename, pedigrees_filename, delta, max_marker_missing_data);

  // *****  done processing command line  *****
  GenotypesSet* the_genotypes_set;
  if(0){
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  long markerid_count = 0;
  Vstr* marker_ids = construct_vstr(1000);
  char* saveptr;
  if((nread = getline(&line, &len, g_stream)) != -1){
    char* token = strtok_r(line, "\t \n", &saveptr);
    // printf("0line (pointer): %p \n", line);
    if((token == NULL)  || (strcmp(token, "MARKER") != 0)){
      exit(EXIT_FAILURE);
    }
    // printf("%s  ", token);
    while(1){
      token = strtok_r(NULL, "\t \n", &saveptr);
      //  printf("line (pointer): %p \n", line);
      if(token == NULL) break;
      char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char)); 
      add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
      //  printf("%s ", token);
      markerid_count++;
    }
    // printf("\n");
  }
  //  free(line); line = NULL;
  fprintf(stderr, "# number of marker ids counted: %ld \n", markerid_count);

  long accession_count = 0;
  Vstr* accession_ids = construct_vstr(1000);
  Vstr* genotypes_strings = construct_vstr(1000);
  long* missing_data_counts = (long*)calloc(markerid_count, sizeof(long));
  Vlong* the_md_vlong = construct_vlong_from_array(markerid_count, missing_data_counts);
  while((nread = getline(&line, &len, g_stream)) != -1){
    accession_count++;
    char* token = strtok_r(line, "\t \n", &saveptr);
    //printf("accession number: %ld id: %s\n", accession_count, token);
    char* acc_id = (char*)malloc((strlen(token)+1)*sizeof(char)); 
    add_string_to_vstr(accession_ids, strcpy(acc_id, token)); // store copy of accession id
    //  printf("%s  ", token);
    long marker_count = 0;
    char* genotypes = (char*)malloc((markerid_count+1) * sizeof(char));
    genotypes[markerid_count] = '\0'; // terminate with null.
    while(1){
      token = strtok_r(NULL, "\t \n", &saveptr);
      if(token == NULL) break;
      double dosage = atof(token);

      if(dosage <= delta){
	genotypes[marker_count] = '0';
      }else if((dosage >= 1-delta) && (dosage <= 1+delta)){
	genotypes[marker_count] = '1';
      }else if(dosage >= 2-delta){
	genotypes[marker_count] = '2';
      }else{ // missing data
	genotypes[marker_count] = '3';
	missing_data_counts[marker_count]++;
      }
      marker_count++;
    } // done reading dosages for all markers, and rounding to 0, 1, 2, 3
    add_string_to_vstr(genotypes_strings, genotypes);
    if(marker_count != markerid_count) exit(EXIT_FAILURE);
    // printf("%s\n", genotypes);
  } // done reading all lines  
  free(line); // only needs to be freed once.
  fclose(g_stream);
  
  the_genotypes_set = construct_genotypesset(accession_ids, marker_ids, genotypes_strings, the_md_vlong);
  }else{
    the_genotypes_set = read_genotypes_file_and_store(g_stream, delta, max_marker_missing_data);
    fclose(g_stream);
  }
  // ***************  read the pedigree file  ********************************
  fclose(p_stream);

  // ***************  Done reading input files  ******************************


  //print_genotypesset(the_genotypes_set);
  print_genotypesset_summary_info(the_genotypes_set);
  if(DBUG) check_genotypesset(the_genotypes_set, max_marker_missing_data);
  
  long good_marker_count = 0;
  for(long i=0; i<the_genotypes_set->marker_ids->size; i++){
    // printf("i: %ld  missing data count: %ld \n", i, missing_data_count[i]);
    if((double)(the_genotypes_set->marker_missing_data_counts->a[i])/(double)the_genotypes_set->n_accessions < max_marker_missing_data) good_marker_count++;
  }
  fprintf(stderr, "# number of good markers: %ld\n", good_marker_count);

  GenotypesSet* the_cleaned_genotypes_set = construct_cleaned_genotypesset(the_genotypes_set, max_marker_missing_data);
  free_genotypesset(the_genotypes_set);
    print_genotypesset_summary_info(the_cleaned_genotypes_set);
  print_genotypesset(the_cleaned_genotypes_set);
   if(DBUG) check_genotypesset(the_cleaned_genotypes_set, max_marker_missing_data);



// ********************  cleanup  **************************
   free_genotypesset(the_cleaned_genotypes_set);
   // getchar();
} // ******************  end of main  ***********************


// **********************************************************
// ******************  functions  ***************************
// **********************************************************

GenotypesSet* read_genotypes_file_and_store(FILE* g_stream, double delta, double max_missing_data_fraction){
 char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  long markerid_count = 0;
  Vstr* marker_ids = construct_vstr(1000);
  char* saveptr;
  if((nread = getline(&line, &len, g_stream)) != -1){
    char* token = strtok_r(line, "\t \n", &saveptr);
    // printf("0line (pointer): %p \n", line);
    if((token == NULL)  || (strcmp(token, "MARKER") != 0)){
      exit(EXIT_FAILURE);
    }
    // printf("%s  ", token);
    while(1){
      token = strtok_r(NULL, "\t \n", &saveptr);
      //  printf("line (pointer): %p \n", line);
      if(token == NULL) break;
      char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char)); 
      add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
      //  printf("%s ", token);
      markerid_count++;
    }
    // printf("\n");
  }
  //  free(line); line = NULL;
  fprintf(stderr, "# number of marker ids counted: %ld \n", markerid_count);

  long accession_count = 0;
  Vstr* accession_ids = construct_vstr(1000);
  Vstr* genotypes_strings = construct_vstr(1000);
  long* missing_data_counts = (long*)calloc(markerid_count, sizeof(long));
  Vlong* the_md_vlong = construct_vlong_from_array(markerid_count, missing_data_counts);
  while((nread = getline(&line, &len, g_stream)) != -1){
    accession_count++;
    char* token = strtok_r(line, "\t \n", &saveptr);
    //printf("accession number: %ld id: %s\n", accession_count, token);
    char* acc_id = (char*)malloc((strlen(token)+1)*sizeof(char)); 
    add_string_to_vstr(accession_ids, strcpy(acc_id, token)); // store copy of accession id
    //  printf("%s  ", token);
    long marker_count = 0;
    char* genotypes = (char*)malloc((markerid_count+1) * sizeof(char));
    genotypes[markerid_count] = '\0'; // terminate with null.
    while(1){
      token = strtok_r(NULL, "\t \n", &saveptr);
      if(token == NULL) break;
      double dosage = atof(token);

      if(dosage <= delta){
	genotypes[marker_count] = '0';
      }else if((dosage >= 1-delta) && (dosage <= 1+delta)){
	genotypes[marker_count] = '1';
      }else if(dosage >= 2-delta){
	genotypes[marker_count] = '2';
      }else{ // missing data
	genotypes[marker_count] = '3';
	missing_data_counts[marker_count]++;
      }
      marker_count++;
    } // done reading dosages for all markers, and rounding to 0, 1, 2, 3
    add_string_to_vstr(genotypes_strings, genotypes);
    if(marker_count != markerid_count) exit(EXIT_FAILURE);
    // printf("%s\n", genotypes);
  } // done reading all lines
  
  free(line); // only needs to be freed once.
  GenotypesSet* the_genotypes_set = construct_genotypesset(accession_ids, marker_ids, genotypes_strings, the_md_vlong);
  return the_genotypes_set;
}
